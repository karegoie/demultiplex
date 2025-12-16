use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;

use anyhow::{anyhow, Context, Result};
use bio::alphabets::dna::revcomp;
use clap::Parser;
use flate2::read::MultiGzDecoder;
use rayon::prelude::*;
use seq_io::fastq::{Reader, Record};
use serde::Serialize;

#[derive(Parser, Debug)]
#[command(author, version, about = "Multiplexed amplicon demultiplexing and quantification")]
struct Args {
    /// Path to R1 FASTQ(.gz)
    #[arg(long)]
    r1: PathBuf,
    /// Path to R2 FASTQ(.gz)
    #[arg(long)]
    r2: PathBuf,
    /// CSV with Sample_ID,Primer_A,Primer_B
    #[arg(long)]
    samples: PathBuf,
    /// CSV/FASTA with Gene_ID,WT_Sequence
    #[arg(long)]
    genes: PathBuf,
    /// Output CSV path
    #[arg(long)]
    output: PathBuf,
    /// Number of threads
    #[arg(long, default_value_t = num_cpus::get())]
    threads: usize,
}

#[derive(Debug, Clone)]
struct SampleRecord {
    id: String,
    primer_a: String,
    primer_b: String,
    primer_a_rc: String,
    primer_b_rc: String,
}

#[derive(Debug, Clone)]
struct GeneRecord {
    id: String,
    wt: String,
}

#[derive(Serialize)]
struct OutputRow {
    sample_id: String,
    gene_id: String,
    r#type: String,
    sequence: String,
    count: u64,
    percentage: f64,
}

#[derive(Copy, Clone, Debug)]
enum Orientation {
    Forward,
    Swapped,
}

fn read_samples(path: &PathBuf) -> Result<Vec<SampleRecord>> {
    let mut rdr = csv::Reader::from_path(path)?;
    let mut out = Vec::new();
    for rec in rdr.records() {
        let rec = rec?;
        let sample_id = rec
            .get(0)
            .ok_or_else(|| anyhow!("Sample_ID missing"))?
            .to_string();
        let primer_a = rec
            .get(1)
            .ok_or_else(|| anyhow!("Primer_A missing"))?
            .to_ascii_uppercase();
        let primer_b = rec
            .get(2)
            .ok_or_else(|| anyhow!("Primer_B missing"))?
            .to_ascii_uppercase();
        out.push(SampleRecord {
            id: sample_id,
            primer_a_rc: String::from_utf8(revcomp(primer_a.as_bytes()))?,
            primer_b_rc: String::from_utf8(revcomp(primer_b.as_bytes()))?,
            primer_a,
            primer_b,
        });
    }
    Ok(out)
}

fn read_genes(path: &PathBuf) -> Result<Vec<GeneRecord>> {
    // CSV reader (Gene_ID, WT_Sequence). Extend to FASTA if desired.
    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(true)
        .flexible(true)
        .from_path(path)?;
    let mut genes = Vec::new();
    for rec in rdr.records() {
        let rec = rec?;
        let id = rec
            .get(0)
            .ok_or_else(|| anyhow!("Gene_ID missing"))?
            .to_string();
        let seq = rec
            .get(1)
            .ok_or_else(|| anyhow!("WT_Sequence missing"))?
            .to_ascii_uppercase();
        genes.push(GeneRecord { id, wt: seq });
    }
    Ok(genes)
}

fn load_fastq_pairs(r1: &PathBuf, r2: &PathBuf) -> Result<Vec<(Vec<u8>, Vec<u8>)>> {
    let r1_file = File::open(r1).with_context(|| format!("Opening {:?}", r1))?;
    let r2_file = File::open(r2).with_context(|| format!("Opening {:?}", r2))?;
    let mut r1_reader = Reader::new(MultiGzDecoder::new(BufReader::new(r1_file)));
    let mut r2_reader = Reader::new(MultiGzDecoder::new(BufReader::new(r2_file)));

    let mut pairs = Vec::new();
    loop {
        let r1_rec = match r1_reader.next() {
            Some(res) => res?,
            None => break,
        };
        let r2_rec = match r2_reader.next() {
            Some(res) => res?,
            None => break,
        };
        pairs.push((r1_rec.seq().to_vec(), r2_rec.seq().to_vec()));
    }
    Ok(pairs)
}

fn detect_sample<'a>(
    samples: &'a [SampleRecord],
    r1: &[u8],
    r2: &[u8],
) -> Option<(&'a SampleRecord, Orientation)> {
    for s in samples {
        if r1.starts_with(s.primer_a.as_bytes()) && r2.starts_with(s.primer_b.as_bytes()) {
            return Some((s, Orientation::Forward));
        }
        if (r1.starts_with(s.primer_b.as_bytes()) || r1.starts_with(s.primer_b_rc.as_bytes()))
            && (r2.starts_with(s.primer_a.as_bytes()) || r2.starts_with(s.primer_a_rc.as_bytes()))
        {
            return Some((s, Orientation::Swapped));
        }
    }
    None
}

fn hamming(a: &[u8], b: &[u8]) -> Option<usize> {
    if a.len() != b.len() {
        return None;
    }
    Some(a.iter().zip(b.iter()).filter(|(x, y)| x != y).count())
}

fn assign_gene<'a>(genes: &'a [GeneRecord], seq: &[u8]) -> (&'a str, bool) {
    let mut best: Option<(&GeneRecord, usize)> = None;
    for g in genes {
        if let Some(dist) = hamming(seq, g.wt.as_bytes()) {
            match best {
                None => best = Some((g, dist)),
                Some((_, d_old)) if dist < d_old => best = Some((g, dist)),
                _ => {}
            }
        }
    }
    if let Some((g, dist)) = best {
        (g.id.as_str(), dist == 0)
    } else {
        ("Unknown", false)
    }
}

fn merge_reads(r1: &[u8], r2: &[u8], min_overlap: usize, max_mismatch: usize) -> Option<Vec<u8>> {
    let max_olap = usize::min(r1.len(), r2.len());
    let mut best: Option<(usize, usize)> = None; // (overlap, mismatches)
    for olap in (min_overlap..=max_olap).rev() {
        let mut mismatches = 0;
        for i in 0..olap {
            if r1[r1.len() - olap + i] != r2[i] {
                mismatches += 1;
                if mismatches > max_mismatch {
                    break;
                }
            }
        }
        if mismatches <= max_mismatch {
            best = Some((olap, mismatches));
            break;
        }
    }
    if let Some((olap, _)) = best {
        let mut merged = Vec::with_capacity(r1.len() + r2.len() - olap);
        merged.extend_from_slice(r1);
        merged.extend_from_slice(&r2[olap..]);
        Some(merged)
    } else {
        None
    }
}

fn main() -> Result<()> {
    let args = Args::parse();
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .context("Configuring rayon")?;

    let samples = read_samples(&args.samples).context("Reading samples")?;
    let genes = read_genes(&args.genes).context("Reading genes")?;
    let pairs = load_fastq_pairs(&args.r1, &args.r2).context("Loading FASTQ pairs")?;

    // Process in parallel
    let counts: HashMap<(String, String, String, String), u64> = pairs
        .par_iter()
        .filter_map(|(r1, r2)| {
            let (sample, orientation) = detect_sample(&samples, r1, r2)?;
            let merged = merge_reads(r1, r2, 20, 3)?;
            let oriented = match orientation {
                Orientation::Forward => merged,
                Orientation::Swapped => revcomp(&merged),
            };
            let (gene_id, is_wt) = assign_gene(&genes, &oriented);
            let variant_seq = String::from_utf8_lossy(&oriented).into_owned();
            let variant_type = if is_wt { "WT" } else { "Mutant" }.to_string();
            Some((
                (
                    sample.id.clone(),
                    gene_id.to_string(),
                    variant_type,
                    variant_seq,
                ),
                1u64,
            ))
        })
        .fold(HashMap::new, |mut acc, (key, c)| {
            *acc.entry(key).or_insert(0) += c;
            acc
        })
        .reduce(HashMap::new, |mut a, b| {
            for (k, v) in b {
                *a.entry(k).or_insert(0) += v;
            }
            a
        });

    // Compute totals per Sample_ID + Gene_ID for percentages
    let mut totals: HashMap<(String, String), u64> = HashMap::new();
    for ((sample, gene, _, _), c) in counts.iter() {
        *totals.entry((sample.clone(), gene.clone())).or_insert(0) += *c;
    }

    let mut wtr = csv::Writer::from_path(&args.output)?;
    for ((sample, gene, vtype, seq), count) in counts.iter() {
        let total = *totals.get(&(sample.clone(), gene.clone())).unwrap_or(&0);
        let denom = if total == 0 { 1 } else { total };
        let pct = (*count as f64) * 100.0 / (denom as f64);
        wtr.serialize(OutputRow {
            sample_id: sample.clone(),
            gene_id: gene.clone(),
            r#type: vtype.clone(),
            sequence: seq.clone(),
            count: *count,
            percentage: pct,
        })?;
    }
    wtr.flush()?;
    Ok(())
}
