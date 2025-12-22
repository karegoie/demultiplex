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

// --- 1. Define Reference Data ---

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

#[derive(Copy, Clone, Debug)]
enum ReadOrientation {
    R1FwdR2Rev,
    R1RevR2Fwd,
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
    let r2_rc = bio::alphabets::dna::revcomp(r2); // Use bio crate for revcomp
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

fn calculate_similarity_bytes(seq1: &[u8], seq2: &[u8]) -> f64 {
    let mut matches = 0;
    let min_len = std::cmp::min(seq1.len(), seq2.len());
    if min_len == 0 { return 0.0; }
    for i in 0..min_len {
        if seq1[i] == seq2[i] {
            matches += 1;
        }
    }
    matches as f64 / min_len as f64
}

fn separate_gene_read(sequence: &[u8], genes: &[GeneRecord]) -> Option<(String, Vec<u8>)> {
    const SIMILARITY_THRESHOLD: f64 = 0.8; 

    let mut best_gene_id: Option<&str> = None;
    let mut best_similarity: f64 = -1.0;

    for gene in genes {
        let sim = calculate_similarity_bytes(sequence, gene.wt.as_bytes());
        if sim > best_similarity {
            best_similarity = sim;
            best_gene_id = Some(&gene.id);
        }
    }

    if let Some(gene_id) = best_gene_id {
        if best_similarity >= SIMILARITY_THRESHOLD {
            return Some((gene_id.to_string(), sequence.to_vec()));
        }
    }
    None
}

fn classify_mutation_bytes(sequence: &[u8], wt_sequence: &[u8]) -> (String, Vec<u8>) {
    if sequence == wt_sequence {
        ("Wild Type".to_string(), wt_sequence.to_vec())
    } else {
        ("Mutation".to_string(), sequence.to_vec())
    }
}

// --- 5. Main Pipeline & Output ---

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

const CHUNK_SIZE: usize = 10000; // Process 10000 paired-end reads at a time

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
            (None, None) => break, // End of both files
            (Some(Err(e)), _) | (_, Some(Err(e))) => return Err(e.into()), // Error reading
            (Some(_), None) | (None, Some(_)) => return Err(anyhow!("Mismatched number of reads in R1 and R2 files").into()), // Mismatched reads
        }
    }
    Ok(chunk)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();

    rayon::ThreadPoolBuilder::new().num_threads(num_cpus::get()).build_global()?;

    // --- Read reference data ---
    let samples = read_samples(&args.samples).context("Failed to read samples file")?;
    let genes = read_genes(&args.genes).context("Failed to read genes file")?;

    let gene_map: HashMap<String, String> = genes.iter().map(|g| (g.id.clone(), g.wt.clone())).collect();

    let r1_file = File::open(&args.r1)?;
    let r2_file = File::open(&args.r2)?;

    let mut r1_reader = FastqReader::new(MultiGzDecoder::new(BufReader::new(r1_file)));
    let mut r2_reader = FastqReader::new(MultiGzDecoder::new(BufReader::new(r2_file)));

    let mut results_map: HashMap<(String, String, String, String), u32> = HashMap::new();

    let mut total_reads_processed = 0;
    // Process reads in a streaming and parallel manner
    loop {
        let chunk = read_pairs_chunked(&mut r1_reader, &mut r2_reader)?;
        if chunk.is_empty() { break; }

        let chunk_results: Vec<HashMap<(String, String, String, String), u32>> = chunk.par_iter().map(|(r1_seq, r2_seq)| {
            let mut local_results_map: HashMap<(String, String, String, String), u32> = HashMap::new();
            
            // 1. Demultiplexing and Read Orientation Detection
            if let Some((sample_record, orientation)) = detect_sample(&samples, r1_seq, r2_seq) {
                // 2. Merge Reads
                if let Some(merged_seq) = merge_reads(r1_seq, r2_seq, 10, 3) {
                    // 3. Trimming
                    let mut trimmed_seq = merged_seq.clone();
                    let start_trim = match orientation {
                        ReadOrientation::R1FwdR2Rev => sample_record.primer_a_len,
                        ReadOrientation::R1RevR2Fwd => sample_record.primer_b_len,
                    };
                    let end_trim = match orientation {
                        ReadOrientation::R1FwdR2Rev => sample_record.primer_b_len,
                        ReadOrientation::R1RevR2Fwd => sample_record.primer_a_len,
                    };

                    if trimmed_seq.len() > start_trim + end_trim {
                        trimmed_seq = trimmed_seq[start_trim..trimmed_seq.len() - end_trim].to_vec();
                    } else {
                        // If trimmed sequence is too short, skip
                        return local_results_map; 
                    }
                    
                    // 4. Gene Separation
                    if let Some((gene_id, gene_sequence_bytes)) = separate_gene_read(&trimmed_seq, &genes) {
                        let wt_sequence_str = gene_map.get(&gene_id).unwrap(); 
                        
                        // 5. Mutation Classification
                        let (variant_type_str, sequence_bytes) = classify_mutation_bytes(&gene_sequence_bytes, wt_sequence_str.as_bytes());

                        *local_results_map.entry((sample_record.id.clone(), gene_id, variant_type_str, String::from_utf8_lossy(&sequence_bytes).to_string())).or_insert(0) += 1;
                    } else {
                        // Unassigned gene within a demultiplexed sample
                        *local_results_map.entry((sample_record.id.clone(), "N/A".to_string(), "Unassigned Gene".to_string(), String::from_utf8_lossy(&merged_seq).to_string())).or_insert(0) += 1;
                    }

                } else {
                    // Unmerged read in a demultiplexed sample
                    *local_results_map.entry((sample_record.id.clone(), "N/A".to_string(), "Unmerged".to_string(), format!("R1:{} R2:{}", String::from_utf8_lossy(r1_seq), String::from_utf8_lossy(r2_seq)))).or_insert(0) += 1;
                }
            } else {
                // Unidentified sample
                *local_results_map.entry(("Unidentified".to_string(), "N/A".to_string(), "N/A".to_string(), format!("R1:{} R2:{}", String::from_utf8_lossy(r1_seq), String::from_utf8_lossy(r2_seq)))).or_insert(0) += 1;
            }
            local_results_map
        }).collect();

        // Aggregate local results into the main results_map
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

    // --- Generate Output CSV ---
    let mut wtr = csv::Writer::from_path(&args.output)?;
    for ((sample_id, gene, variant_type, sequence), count) in results_map {
        wtr.serialize(OutputRow {
            sample_id,
            gene,
            variant_type,
            sequence,
            count,
        })?;
    }
    wtr.flush()?;
    println!("Analysis complete. Summary saved to {}", args.output.display());

    Ok(())
}