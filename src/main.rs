use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, Write};
use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering};

use anyhow::{anyhow, Context, Result};
use bio::alphabets::dna::revcomp;
use block_aligner::scan_block::{Block, PaddedBytes};
use block_aligner::scores::{Gaps, NucMatrix};
use clap::Parser;
use flate2::read::MultiGzDecoder;
use rayon::prelude::*;
use seq_io::fastq::{Reader, Record};
use serde::Serialize;

#[derive(Parser, Debug)]
#[command(author, version, about = "Multiplexed amplicon demultiplexing and quantification")]
struct Args {
    #[arg(long)]
    r1: PathBuf,
    #[arg(long)]
    r2: PathBuf,
    #[arg(long)]
    samples: PathBuf,
    #[arg(long)]
    genes: PathBuf,
    #[arg(long)]
    output: PathBuf,
    #[arg(long, default_value_t = num_cpus::get())]
    threads: usize,
    #[arg(long, default_value_t = false)]
    debug: bool,
}

#[derive(Debug, Clone)]
struct SampleRecord {
    id: String,
    primer_a_anchor: Vec<u8>,
    primer_b_anchor: Vec<u8>,
    primer_a_len: usize,
    primer_b_len: usize,
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
enum ReadOrientation {
    R1FwdR2Rev,
    R1RevR2Fwd,
}

fn read_genes(path: &PathBuf) -> Result<Vec<GeneRecord>> {
    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(true)
        .flexible(true)
        .from_path(path)?;
    let mut genes = Vec::new();
    for rec in rdr.records() {
        let rec = rec?;
        let id = rec.get(0).ok_or_else(|| anyhow!("Gene_ID missing"))?.to_string();
        let seq = rec.get(1).ok_or_else(|| anyhow!("WT_Sequence missing"))?.to_ascii_uppercase();
        genes.push(GeneRecord { id, wt: seq });
    }
    Ok(genes)
}

fn find_lcp(strings: &[String]) -> String {
    if strings.is_empty() { return String::new(); }
    let first = &strings[0];
    let mut len = first.len();
    for s in &strings[1..] {
        while len > 0 {
            if s.starts_with(&first[..len]) { break; }
            len -= 1;
        }
    }
    first[..len].to_string()
}

fn read_samples(path: &PathBuf) -> Result<Vec<SampleRecord>> {
    let mut rdr = csv::Reader::from_path(path)?;
    let mut raw_records = Vec::new();
    for rec in rdr.records() {
        let rec = rec?;
        let id = rec.get(0).unwrap().to_string();
        let pa = rec.get(1).unwrap().to_ascii_uppercase();
        let pb = rec.get(2).unwrap().to_ascii_uppercase();
        raw_records.push((id, pa, pb));
    }
    if raw_records.is_empty() { return Ok(Vec::new()); }

    let pa_list: Vec<String> = raw_records.iter().map(|r| r.1.clone()).collect();
    let pb_list: Vec<String> = raw_records.iter().map(|r| r.2.clone()).collect();
    let adapter_a = find_lcp(&pa_list);
    let adapter_b = find_lcp(&pb_list);

    let mut out = Vec::new();
    const ANCHOR_LEN: usize = 16; 
    for (id, pa, pb) in raw_records {
        let target_a = &pa[adapter_a.len()..];
        let target_b = &pb[adapter_b.len()..];
        let anchor_a = if target_a.len() >= ANCHOR_LEN { &target_a[..ANCHOR_LEN] } else { target_a };
        let anchor_b = if target_b.len() >= ANCHOR_LEN { &target_b[..ANCHOR_LEN] } else { target_b };
        out.push(SampleRecord {
            id,
            primer_a_anchor: anchor_a.as_bytes().to_vec(),
            primer_b_anchor: anchor_b.as_bytes().to_vec(),
            primer_a_len: target_a.len(),
            primer_b_len: target_b.len(),
        });
    }
    Ok(out)
}

fn load_fastq_pairs(r1: &PathBuf, r2: &PathBuf) -> Result<Vec<(Vec<u8>, Vec<u8>)>> {
    let r1_file = File::open(r1)?;
    let r2_file = File::open(r2)?;
    let mut r1_reader = Reader::new(MultiGzDecoder::new(BufReader::new(r1_file)));
    let mut r2_reader = Reader::new(MultiGzDecoder::new(BufReader::new(r2_file)));
    let mut pairs = Vec::new();
    while let (Some(r1), Some(r2)) = (r1_reader.next(), r2_reader.next()) {
        pairs.push((r1?.seq().to_vec(), r2?.seq().to_vec()));
    }
    Ok(pairs)
}

fn matches_start(read: &[u8], primer: &[u8], max_mismatch: usize) -> bool {
    if read.len() < primer.len() { return false; }
    read.iter().zip(primer.iter()).filter(|(a, b)| a != b).count() <= max_mismatch
}

fn detect_sample<'a>(samples: &'a [SampleRecord], r1: &[u8], r2: &[u8]) -> Option<(&'a SampleRecord, ReadOrientation)> {
    const MAX_MM: usize = 2; 
    for s in samples {
        if matches_start(r1, &s.primer_a_anchor, MAX_MM) && matches_start(r2, &s.primer_b_anchor, MAX_MM) {
            return Some((s, ReadOrientation::R1FwdR2Rev));
        }
        if matches_start(r1, &s.primer_b_anchor, MAX_MM) && matches_start(r2, &s.primer_a_anchor, MAX_MM) {
            return Some((s, ReadOrientation::R1RevR2Fwd));
        }
    }
    None
}

fn merge_reads(r1: &[u8], r2: &[u8], min_overlap: usize, max_mismatch: usize) -> Option<Vec<u8>> {
    let len1 = r1.len();
    let r2_rc = revcomp(r2);
    let max_olap = len1.min(r2_rc.len());
    for olap in (min_overlap..=max_olap).rev() {
        let r1_suffix = &r1[len1 - olap..];
        let r2_prefix = &r2_rc[..olap];
        if r1_suffix.iter().zip(r2_prefix.iter()).filter(|(a,b)| a!=b).count() <= max_mismatch {
            let mut merged = Vec::with_capacity(len1 + r2_rc.len() - olap);
            merged.extend_from_slice(&r1[..len1-olap]);
            merged.extend_from_slice(&r2_rc);
            return Some(merged);
        }
    }
    None 
}

fn best_alignment_dual<'a>(genes: &'a [GeneRecord], seq: &[u8]) -> (&'a str, bool, f64, String) {
    if genes.is_empty() {
        return ("Unknown", false, f64::NEG_INFINITY, String::from_utf8_lossy(seq).to_string());
    }
    let (id_f, score_f) = align_to_genes(genes, seq);
    let seq_rc = revcomp(seq);
    let (id_r, score_r) = align_to_genes(genes, &seq_rc);

    if score_f >= score_r {
        (id_f, false, score_f, String::from_utf8_lossy(seq).to_string())
    } else {
        (id_r, false, score_r, String::from_utf8_lossy(&seq_rc).to_string())
    }
}

fn align_to_genes<'a>(genes: &'a [GeneRecord], seq: &[u8]) -> (&'a str, f64) {
    let max_ref = genes.iter().map(|g| g.wt.len()).max().unwrap_or(1).max(1);
    let max_query = seq.len().max(1);
    let max_len = max_ref.max(max_query);
    let max_block_size = max_len.next_power_of_two().max(64).min(2048); 
    
    let matrix = NucMatrix::new_simple(2, -3); 
    let gaps = Gaps { open: -5, extend: -1 }; 
    let mut block = Block::<true, false>::new(max_query, max_ref, max_block_size);
    let q = PaddedBytes::from_bytes::<NucMatrix>(seq, max_block_size);

    let mut best_id = "Unknown";
    let mut best_norm = f64::NEG_INFINITY;

    for g in genes {
        let r = PaddedBytes::from_bytes::<NucMatrix>(g.wt.as_bytes(), max_block_size);
        block.align(&q, &r, &matrix, gaps, 32..=max_block_size, 0); 
        let res = block.res();
        let score = res.score as f64;
        let max_score = (2 * g.wt.len().max(seq.len())) as f64; 
        let norm = score / max_score;

        if norm > best_norm {
            best_norm = norm;
            best_id = g.id.as_str();
        }
    }
    (best_id, best_norm)
}

fn main() -> Result<()> {
    let args = Args::parse();
    rayon::ThreadPoolBuilder::new().num_threads(args.threads).build_global()?;

    let genes = read_genes(&args.genes).context("Reading genes")?;
    let gene_map: HashMap<String, String> = genes.iter().map(|g| (g.id.clone(), g.wt.clone())).collect();
    
    let samples = read_samples(&args.samples).context("Reading samples")?;
    let pairs = load_fastq_pairs(&args.r1, &args.r2).context("Loading FASTQ pairs")?;

    let demux_hits = AtomicUsize::new(0);
    let merged_success = AtomicUsize::new(0);

    let debug_log = if args.debug {
        Some(File::create("debug_mismatches.txt")?)
    } else {
        None
    };

    let counts: HashMap<(String, String, String), (String, u64)> = pairs
        .par_iter()
        .filter_map(|(r1, r2)| {
            let (sample, orientation) = detect_sample(&samples, r1, r2)?;
            demux_hits.fetch_add(1, Ordering::Relaxed);

            let merged = merge_reads(r1, r2, 10, 3)?; 
            merged_success.fetch_add(1, Ordering::Relaxed);

            let mut trimmed = merged.clone();
            let start_trim = match orientation {
                ReadOrientation::R1FwdR2Rev => sample.primer_a_len,
                ReadOrientation::R1RevR2Fwd => sample.primer_b_len,
            };
            let end_trim = match orientation {
                ReadOrientation::R1FwdR2Rev => sample.primer_b_len,
                ReadOrientation::R1RevR2Fwd => sample.primer_a_len,
            };

            if trimmed.len() > start_trim + end_trim + 10 {
                trimmed = trimmed[start_trim..trimmed.len() - end_trim].to_vec();
            } else {
                return None; 
            }

            let (gene_id, _, _score, oriented_seq) = best_alignment_dual(&genes, &trimmed);
            
            let wt_seq_str = gene_map.get(gene_id).unwrap();
            let seq_upper = oriented_seq.to_ascii_uppercase();
            let wt_upper = wt_seq_str.to_ascii_uppercase();

            let is_contained = wt_upper.contains(&seq_upper) || seq_upper.contains(&wt_upper);
            let final_wt_check = is_contained; 
            
            let type_str = if final_wt_check { "WT" } else { "Mutant" };
            
            let final_seq = if final_wt_check { 
                wt_seq_str.clone()
            } else { 
                oriented_seq.clone() 
            };

            Some(((sample.id.clone(), gene_id.to_string(), final_seq), (type_str.to_string(), 1u64)))
        })
        .fold(HashMap::new, |mut acc, (key, (vtype, c))| {
            let entry = acc.entry(key).or_insert((vtype, 0));
            entry.1 += c;
            acc
        })
        .reduce(HashMap::new, |mut a, b| {
            for (k, (vtype, count)) in b {
                let entry = a.entry(k).or_insert((vtype, 0));
                entry.1 += count;
            }
            a
        });

    let mut wtr = csv::Writer::from_path(&args.output)?;
    
    let mut group_totals: HashMap<(String, String), u64> = HashMap::new();
    for ((s, g, _seq), (_, c)) in counts.iter() {
        *group_totals.entry((s.clone(), g.clone())).or_insert(0) += *c;
    }

    let mut rows: Vec<_> = counts.into_iter().collect();
    rows.sort_by(|a, b| {
        a.0.0.cmp(&b.0.0)
            .then(a.0.1.cmp(&b.0.1))
            .then(b.1.0.cmp(&a.1.0))
            .then(b.1.1.cmp(&a.1.1))
    });

    if let Some(mut f) = debug_log {
        writeln!(f, "SampleID\tGeneID\tCount\tType\tRead_Sequence\tExpected_WT_Sequence")?;
        for ((s, g, seq), (t, c)) in rows.iter().take(50) {
            if t == "Mutant" {
                let expected = gene_map.get(g).unwrap();
                writeln!(f, "{}\t{}\t{}\t{}\t{}\t{}", s, g, c, t, seq, expected)?;
            }
        }
        eprintln!("Debug log written to 'debug_mismatches.txt'.");
    }

    for ((sample, gene, seq), (vtype, count)) in rows {
        let total = *group_totals.get(&(sample.clone(), gene.clone())).unwrap_or(&1);
        let pct = (count as f64 / total as f64) * 100.0;
        
        wtr.serialize(OutputRow {
            sample_id: sample,
            gene_id: gene,
            r#type: vtype,
            sequence: seq,
            count,
            percentage: pct,
        })?;
    }
    wtr.flush()?;

    if args.debug {
        eprintln!("Demultiplexed Reads: {}", demux_hits.load(Ordering::Relaxed));
        eprintln!("Successfully Merged: {}", merged_success.load(Ordering::Relaxed));
    }

    Ok(())
}
