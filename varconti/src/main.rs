#![allow(dead_code)]

use std::time::{Instant};

mod parsefa;
mod splitter;
mod align;
mod cliparameters;

use std::fs::File;
use std::io::{Error, Write};

extern crate bincode;

extern crate chrono;
use chrono::Local;

mod serializer;
mod indexscan;

use colored::*;

use std::hash::{Hash, Hasher};
use std::collections::HashMap;
use std::collections::hash_map::DefaultHasher;


// cargo build --release
// "target/release/varconti.exe" --index -f "files/Homo_sapiens.GRCh38.cdna.all.fa" -k 22 -o "output/myindex_24.idx"
// "target/release/varconti.exe" -f "files/SRR3534129.fastq" -o "output/counts2.tsv -x "output/myindex_24.idx"
// "target/release/varconti.exe" -f "files/SRR3534129.fastq" -o "output/counts2.tsv" -x "output/myindex_24.idx"

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

    let matches = cliparameters::cli();

    // IO parameters
    let input_file = matches.value_of("input-file").unwrap();
    let output_file = matches.value_of("output-file").unwrap();
    let index_file = matches.value_of("index-file").unwrap_or("none");

    let thread_count = matches.value_of("thread-number").unwrap_or("1").parse::<usize>().unwrap();
    let kmer_length = matches.value_of("k-mer").unwrap_or("22").parse::<u32>().unwrap();
    let verbose = matches.is_present("verbose");
    let build_index = matches.is_present("index");

    let now = Instant::now();

    //let kmer_length: u32 = 40;
    let kmer_offset: u32 = 0;
    let step_size = 1;
    let version: u32 = 1;
    let sensitivity = 0;

    let genome_file = "files/Homo_sapiens.GRCh38.cdna.all.fa";
    //let output_idx = "output/myindex_24.idx";
    let output_counts = "output/counts.tsv";
    //let fastq_file = "files/SRR3534129.fastq";
    //let fastq_file = "files/SRR1525256.fastq";
    //let build_index = true;

    if build_index{
        println!("{}", format!("[{}] Build index | {} | k={}", Local::now().format("%Y-%m-%d][%H:%M:%S"), input_file, kmer_length).blue());
        
        let (transcripts, eq_classes, eq_elements, transcript_kmers) = parsefa::read_fa(input_file, kmer_length);

        println!("{}", format!("[{}] Save to file | {}", Local::now().format("%Y-%m-%d][%H:%M:%S"), output_file).blue());
        serializer::serialize(&(output_file.to_string()), &transcripts, &eq_classes, &eq_elements, &transcript_kmers, &kmer_length, &version)?;
        println!("{}", format!("[{}] Done! Elapsed time: {}m {}s", Local::now().format("%Y-%m-%d][%H:%M:%S"), now.elapsed().as_millis()/60000, (now.elapsed().as_millis()%60000)/1000).green());
        println!("{}", format!("Transcripts: {} | EQ Classes: {} | EQ Elements: {}", transcripts.len(), eq_classes.len(), eq_elements.len()).green());
        
        std::thread::spawn(move || drop(eq_classes));
        std::thread::spawn(move || drop(eq_elements));
    }
    else {
        let fastq_file = input_file;
        let output_counts = output_file;
        let now = Instant::now();
        
        println!("{}", format!("[{}] Load index | {}", Local::now().format("%Y-%m-%d][%H:%M:%S"), index_file).blue());
        let (kmer_length, transcripts, eq_elements, eq_classes, transcript_kmers) = serializer::deserialize(&(index_file.to_string()));
        println!("[{}] Index | k: {} | T: {} | EQC: {} | EQE: {}", Local::now().format("%Y-%m-%d][%H:%M:%S"), kmer_length, transcripts.len(), eq_classes.len(), eq_elements.len());
        
        indexscan::index_stats(&transcripts, &eq_classes, &eq_elements, &transcript_kmers, &kmer_length);

        println!("{}", format!("[{}] Align reads | {} | k={}", Local::now().format("%Y-%m-%d][%H:%M:%S"), fastq_file, kmer_length).blue());
        let (line_count, transcript_counts_unique, transcript_counts, unique_count) = align::read_fastq(fastq_file, kmer_length, &eq_classes, &eq_elements, &kmer_offset, step_size, sensitivity);
        
        std::thread::spawn(move || drop(eq_classes));
        std::thread::spawn(move || drop(eq_elements));

        println!("tcl: {}", transcript_counts.len());

        let mut total_count = 0;
        let mut output = File::create(output_counts)?;
        for i in 0..transcripts.len() {
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
    }
    Ok(())
}

