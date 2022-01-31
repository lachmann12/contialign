#![allow(dead_code)]

use std::time::{Instant};

mod parsefa;
mod splitter;
mod align;

use std::fs::File;
use std::io::{Error, Write};

extern crate bincode;

extern crate chrono;
use chrono::Local;

mod serializer;
use colored::*;

use std::hash::{Hash, Hasher};
use std::collections::HashMap;
use std::collections::hash_map::DefaultHasher;


// cargo build --release
// "target/release/varconti.exe"
// "external/kallisto.exe" index -i "output/kallisto/human.idx" "files/Homo_sapiens.GRCh38.cdna.all.fa"
// "external/kallisto.exe" quant -i "output/kallisto/human.idx" -o "output/kallisto/" --single -l 200 -s 20 "files/SRR3534129.fastq"
// "external/kallisto.exe" quant -i "output/kallisto/human.idx" -o "output/kallisto/" --single -l 200 -s 20 "files/SRR1525256.fastq"

pub fn listhash(input: &Vec<u32>) -> u32 {
    let mut h = DefaultHasher::new();
    Hash::hash_slice(input, &mut h);
    return h.finish() as u32;
}

pub fn hash(input: &str) -> u32 {
    let mut hasher = DefaultHasher::new();
    input.hash(&mut hasher);
    return hasher.finish() as u32;
}

fn main() -> Result<(), Error> {

    let now = Instant::now();

    let kmer_length: u32 = 22;
    let kmer_offset: u32 = 0;
    let step_size = 1;
    let version: u32 = 1;
    let sensitivity = 0;

    let genome_file = "files/Homo_sapiens.GRCh38.cdna.all.fa";
    let output_idx = "output/myindex_31.idx";
    let output_counts = "output/counts.tsv";
    //let fastq_file = "files/SRR3534129.fastq";
    let fastq_file = "files/SRR1525256.fastq";
    let build_index = true;

    if build_index{
        println!("{}", format!("[{}] Build index | {} | k={}", Local::now().format("%Y-%m-%d][%H:%M:%S"), genome_file, kmer_length).blue());
        
        let (transcripts, eq_classes, eq_elements) = parsefa::read_fa(genome_file, kmer_length);
        
        println!("{}", format!("[{}] Save to file | {}", Local::now().format("%Y-%m-%d][%H:%M:%S"), output_idx).blue());
        serializer::serialize(&(output_idx.to_string()), &transcripts, &eq_classes, &eq_elements, &kmer_length, &version)?;
        println!("{}", format!("[{}] Done! Elapsed time: {}m {}s", Local::now().format("%Y-%m-%d][%H:%M:%S"), now.elapsed().as_millis()/60000, (now.elapsed().as_millis()%60000)/1000).green());
        println!("{}", format!("Transcripts: {} | EQ Classes: {} | EQ Elements: {}", transcripts.len(), eq_classes.len(), eq_elements.len()).green());
        
        std::thread::spawn(move || drop(eq_classes));
        std::thread::spawn(move || drop(eq_elements));
    }

    let now = Instant::now();
    
    println!("{}", format!("[{}] Load index | {}", Local::now().format("%Y-%m-%d][%H:%M:%S"), output_idx).blue());
    let (kmer_length, transcripts, eq_elements, eq_classes) = serializer::deserialize(&(output_idx.to_string()));
    println!("[{}] Index | k: {} | T: {} | EQC: {} | EQE: {}", Local::now().format("%Y-%m-%d][%H:%M:%S"), kmer_length, transcripts.len(), eq_classes.len(), eq_elements.len());
    
    println!("{}", format!("[{}] Align reads | {} | k={}", Local::now().format("%Y-%m-%d][%H:%M:%S"), fastq_file, kmer_length).blue());
    let (line_count, transcript_counts, unique_count) = align::read_fastq(fastq_file, kmer_length, &eq_classes, &eq_elements, &kmer_offset, step_size, sensitivity);
    
    std::thread::spawn(move || drop(eq_classes));
    std::thread::spawn(move || drop(eq_elements));

    let mut total_count = 0;
    let mut output = File::create(output_counts)?;
    for i in 0..transcripts.len() {
    //for i in 0..10 {
        if transcript_counts.contains_key(&(i as u32)) {
            total_count = total_count + transcript_counts.get(&(i as u32)).unwrap();
            write!(output, "{}\t{}\n", transcripts[i], transcript_counts.get(&(i as u32)).unwrap())?;
        }
        else{
            write!(output, "{}\t{}\n", transcripts[i], 0)?;
        }
    }

    println!("{}", format!("[{}] Done! Elapsed time: {}m {}s", Local::now().format("%Y-%m-%d][%H:%M:%S"), now.elapsed().as_millis()/60000, (now.elapsed().as_millis()%60000)/1000).green());
    println!("{}", format!("Transcripts: {} | Reads: {} | Aligned reads: {} | Unique Reads: {}", transcripts.len(), line_count/4, total_count, unique_count).green());
    Ok(())
}

