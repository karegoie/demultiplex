use std::collections::HashMap;
use serde::Serialize;
use clap::Parser;
use std::path::PathBuf;
use std::fs::File;
use flate2::read::MultiGzDecoder;
use std::io::BufReader;
use anyhow::{anyhow, Context, Result};
use seq_io::fastq::{Reader as FastqReader, Record as FastqRecord};
use rayon::prelude::*;
use bio::alphabets::dna::revcomp;
use bio::alignment::pairwise::Aligner;
use block_aligner::scan_block::{Block, PaddedBytes};
use block_aligner::scores::{Gaps, NucMatrix};


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
    #[arg(long, default_value_t = false)]
    debug: bool,
}

#[derive(Debug, Clone)]
struct SampleRecord {
    id: String,
    primer_a: Vec<u8>,
    primer_b: Vec<u8>,
    primer_a_rc: Vec<u8>,
    primer_b_rc: Vec<u8>,
    primer_a_len: usize,
    primer_b_len: usize,
}

#[derive(Debug, Clone)]
struct GeneRecord {
    id: String,
    wt: String,
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


fn longest_common_substring(strings: &[String]) -> String {
    if strings.is_empty() {
        return String::new();
    }
    let first = &strings[0];
    let mut longest = "";

    for i in 0..first.len() {
        for j in i + 1..=first.len() {
            let sub = &first[i..j];
            if strings.iter().all(|s| s.contains(sub)) {
                if sub.len() > longest.len() {
                    longest = sub;
                }
            }
        }
    }
    longest.to_string()
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
    
    let lcs_a = longest_common_substring(&pa_list);
    let lcs_b = longest_common_substring(&pb_list);
    
    eprintln!("LCS A: {}", lcs_a);
    eprintln!("LCS B: {}", lcs_b);

    let mut out = Vec::new();
    for (id, pa, pb) in raw_records {
        let pa_specific = pa.replace(&lcs_a, "");
        let pb_specific = pb.replace(&lcs_b, "");
        
        eprintln!("Sample {}: PA_Spec: {}, PB_Spec: {}", id, pa_specific, pb_specific);

        let primer_a_bytes = pa_specific.as_bytes().to_vec();
        let primer_b_bytes = pb_specific.as_bytes().to_vec();
        out.push(SampleRecord {
            id,
            primer_a_len: primer_a_bytes.len(),
            primer_b_len: primer_b_bytes.len(),
            primer_a: primer_a_bytes.clone(),
            primer_b: primer_b_bytes.clone(),
            primer_a_rc: revcomp(&primer_a_bytes),
            primer_b_rc: revcomp(&primer_b_bytes),
        });
    }
    Ok(out)
}

#[derive(Copy, Clone, Debug)]
enum ReadOrientation {
    R1FwdR2Rev,
    R1RevR2Fwd,
}

fn align_primer<'a>(read: &'a [u8], primer: &'a [u8]) -> Option<(i32, usize)> {
    let score = |a: u8, b: u8| if a == b { 2i32 } else { -3i32 };
    let mut aligner = Aligner::with_capacity(primer.len(), read.len(), -5, -1, &score);
    let alignment = aligner.semiglobal(primer, read);
    
    Some((alignment.score, alignment.yend))
}

fn detect_sample_with_alignment<'a>(samples: &'a [SampleRecord], r1: &'a [u8], r2: &'a [u8], debug: bool) -> Option<(&'a SampleRecord, ReadOrientation, usize, usize)> {
    const MIN_NORMALIZED_SCORE: f64 = 0.6;

    for s in samples {
        // R1-Fwd, R2-Fwd
        if let (Some((score1, end1)), Some((score2, end2))) = (align_primer(r1, &s.primer_a), align_primer(r2, &s.primer_b)) {
            let norm_score1 = score1 as f64 / (2.0 * s.primer_a.len() as f64);
            let norm_score2 = score2 as f64 / (2.0 * s.primer_b.len() as f64);
            if debug { eprintln!("Sample: {}, R1Fwd/R2Fwd scores: {:.2} ({}), {:.2} ({})", s.id, norm_score1, score1, norm_score2, score2); }
            if norm_score1 >= MIN_NORMALIZED_SCORE && norm_score2 >= MIN_NORMALIZED_SCORE {
                return Some((s, ReadOrientation::R1FwdR2Rev, end1, end2));
            }
        }
        // R1-Rev, R2-Rev
        if let (Some((score1, end1)), Some((score2, end2))) = (align_primer(r1, &s.primer_b), align_primer(r2, &s.primer_a)) {
            let norm_score1 = score1 as f64 / (2.0 * s.primer_b.len() as f64);
            let norm_score2 = score2 as f64 / (2.0 * s.primer_a.len() as f64);
            if debug { eprintln!("Sample: {}, R1Rev/R2Rev scores: {:.2}, {:.2}", s.id, norm_score1, norm_score2); }
            if norm_score1 >= MIN_NORMALIZED_SCORE && norm_score2 >= MIN_NORMALIZED_SCORE {
                return Some((s, ReadOrientation::R1RevR2Fwd, end1, end2));
            }
        }
        // R1-Fwd-RC, R2-Fwd-RC
        if let (Some((score1, end1)), Some((score2, end2))) = (align_primer(r1, &s.primer_a_rc), align_primer(r2, &s.primer_b_rc)) {
            let norm_score1 = score1 as f64 / (2.0 * s.primer_a_rc.len() as f64);
            let norm_score2 = score2 as f64 / (2.0 * s.primer_b_rc.len() as f64);
            if debug { eprintln!("Sample: {}, R1FwdRC/R2FwdRC scores: {:.2}, {:.2}", s.id, norm_score1, norm_score2); }
            if norm_score1 >= MIN_NORMALIZED_SCORE && norm_score2 >= MIN_NORMALIZED_SCORE {
                return Some((s, ReadOrientation::R1FwdR2Rev, end1, end2));
            }
        }
        // R1-Rev-RC, R2-Rev-RC
        if let (Some((score1, end1)), Some((score2, end2))) = (align_primer(r1, &s.primer_b_rc), align_primer(r2, &s.primer_a_rc)) {
            let norm_score1 = score1 as f64 / (2.0 * s.primer_b_rc.len() as f64);
            let norm_score2 = score2 as f64 / (2.0 * s.primer_a_rc.len() as f64);
            if debug { eprintln!("Sample: {}, R1RevRC/R2RevRC scores: {:.2}, {:.2}", s.id, norm_score1, norm_score2); }
            if norm_score1 >= MIN_NORMALIZED_SCORE && norm_score2 >= MIN_NORMALIZED_SCORE {
                return Some((s, ReadOrientation::R1RevR2Fwd, end1, end2));
            }
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

fn best_alignment_dual<'a>(genes: &'a [GeneRecord], seq: &[u8]) -> (&'a str, f64, bool) {
    if genes.is_empty() {
        return ("Unknown", f64::NEG_INFINITY, false);
    }
    let (id_f, score_f) = align_to_genes(genes, seq);
    let seq_rc = revcomp(seq);
    let (id_r, score_r) = align_to_genes(genes, &seq_rc);

    if score_f >= score_r {
        (id_f, score_f, false)
    } else {
        (id_r, score_r, true)
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

#[derive(Serialize)]
struct OutputRow {
    #[serde(rename = "Sample_ID")]
    sample_id: String,
    #[serde(rename = "Gene")]
    gene: String,
    #[serde(rename = "Variant_Type")]
    variant_type: String,
    #[serde(rename = "Sequence")]
    sequence: String,
    #[serde(rename = "Count")]
    count: u32,
}

const CHUNK_SIZE: usize = 10000;

fn read_pairs_chunked(
    r1_reader: &mut FastqReader<MultiGzDecoder<BufReader<File>>>,
    r2_reader: &mut FastqReader<MultiGzDecoder<BufReader<File>>>,
) -> Result<Vec<(Vec<u8>, Vec<u8>)>, Box<dyn std::error::Error>> {
    let mut chunk = Vec::with_capacity(CHUNK_SIZE);
    for _ in 0..CHUNK_SIZE {
        match (r1_reader.next(), r2_reader.next()) {
            (Some(Ok(r1_rec)), Some(Ok(r2_rec))) => {
                chunk.push((r1_rec.seq().to_vec(), r2_rec.seq().to_vec()));
            },
            (None, None) => break,
            (Some(Err(e)), _) | (_, Some(Err(e))) => return Err(e.into()),
            (Some(_), None) | (None, Some(_)) => return Err(anyhow!("Mismatched number of reads in R1 and R2 files").into()),
        }
    }
    Ok(chunk)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();
    rayon::ThreadPoolBuilder::new().num_threads(num_cpus::get()).build_global()?;

    let samples = read_samples(&args.samples).context("Failed to read samples file")?;
    let genes = read_genes(&args.genes).context("Failed to read genes file")?;
    let gene_map: HashMap<String, String> = genes.iter().map(|g| (g.id.clone(), g.wt.clone())).collect();

    let r1_file = File::open(&args.r1)?;
    let r2_file = File::open(&args.r2)?;

    let mut r1_reader = FastqReader::new(MultiGzDecoder::new(BufReader::new(r1_file)));
    let mut r2_reader = FastqReader::new(MultiGzDecoder::new(BufReader::new(r2_file)));

    let mut results_map: HashMap<(String, String, String, String), u32> = HashMap::new();
    let mut total_reads_processed = 0;

    loop {
        let chunk = read_pairs_chunked(&mut r1_reader, &mut r2_reader)?;
        if chunk.is_empty() { break; }

        let chunk_results: Vec<HashMap<_, _>> = chunk.par_iter().map(|(r1_seq, r2_seq)| {
            let mut local_results_map = HashMap::new();
            if let Some((sample_record, orientation, end1, end2)) = detect_sample_with_alignment(&samples, r1_seq, r2_seq, args.debug) {
                if let Some(merged_seq) = merge_reads(r1_seq, r2_seq, 10, 5) {
                    let start_trim = end1;
                    let end_trim = end2;
                    if merged_seq.len() > start_trim + end_trim {
                        let trimmed_seq = &merged_seq[start_trim..merged_seq.len() - end_trim];
                        let (gene_id, _score, is_rc) = best_alignment_dual(&genes, trimmed_seq);
                        if gene_id != "Unknown" {
                            let wt_sequence_str = gene_map.get(gene_id).unwrap();
                            let final_gene_seq_bytes = if is_rc { revcomp(trimmed_seq) } else { trimmed_seq.to_vec() };
                            let final_gene_seq = String::from_utf8_lossy(&final_gene_seq_bytes).to_string();
                            let variant_type = if final_gene_seq == *wt_sequence_str { "Wild Type" } else { "Mutation" };
                            *local_results_map.entry((sample_record.id.clone(), gene_id.to_string(), variant_type.to_string(), final_gene_seq)).or_insert(0) += 1;
                        } else {
                            *local_results_map.entry((sample_record.id.clone(), "N/A".to_string(), "Unassigned Gene".to_string(), String::from_utf8_lossy(&merged_seq).to_string())).or_insert(0) += 1;
                        }
                    }
                } else {
                    *local_results_map.entry((sample_record.id.clone(), "N/A".to_string(), "Unmerged".to_string(), format!("R1:{} R2:{}", String::from_utf8_lossy(r1_seq), String::from_utf8_lossy(r2_seq)))).or_insert(0) += 1;
                }
            } else {
                *local_results_map.entry(("Unidentified".to_string(), "N/A".to_string(), "N/A".to_string(), format!("R1:{} R2:{}", String::from_utf8_lossy(r1_seq), String::from_utf8_lossy(r2_seq)))).or_insert(0) += 1;
            }
            local_results_map
        }).collect();

        for local_map in chunk_results {
            for (key, count) in local_map {
                *results_map.entry(key).or_insert(0) += count;
            }
        }

        total_reads_processed += chunk.len();
        if args.debug {
            eprintln!("Processed {} reads...", total_reads_processed);
        }
    }

    let mut wtr = csv::Writer::from_path(&args.output)?;
    for ((sample_id, gene, variant_type, sequence), count) in results_map {
        wtr.serialize(OutputRow { sample_id, gene, variant_type, sequence, count })?;
    }
    wtr.flush()?;
    println!("Analysis complete. Summary saved to {}", args.output.display());

    Ok(())
}