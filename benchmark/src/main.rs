use std::collections::HashMap;
use nohash_hasher::IntMap;
use rustc_hash::FxHashMap;

use std::io::Cursor;

use std::time::{Instant};
use chrono::Local;
use colored::*;
use std::time::Duration;


use std::fs::File;
use std::io::{BufRead, BufReader};
use pbr::ProgressBar;

use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};




use cityhash::cityhash_1_1_1::city_hash_64;

mod cliparameters;

// cargo build --release
// "target/release/benchmark"

pub fn read_fastq(input_file: &str, hm_index: IntMap<u32, u64>, kmer_length: u32) {
    let mut seq: String;
    let mut counter = 0;
    let mut read_length = 0;
    
    let now = Instant::now();

    let mut line_count: i64 = 0;
    let f = BufReader::new(File::open(input_file).unwrap());
    for line in f.lines() {
        line_count = line_count+1;

        if line_count == 2 {
            read_length = line.unwrap().to_string().len();
        }
    }

    let elapsed_read = now.elapsed().as_millis();
    //println!("Lean scan {}", elapsed_read);
    let mut pb = ProgressBar::new(line_count as u64/4);
    pb.set_max_refresh_rate(Some(Duration::from_millis(100)));
    pb.message("     -> aligning | ");

    let mut sum = 0;

    let now = Instant::now();
    let mut total_hashed = 0;
    let f = BufReader::new(File::open(input_file).unwrap());
    for (_i, line) in f.lines().enumerate() {
        
        if counter % 4 == 1 {
            pb.inc();
            seq = line.unwrap().to_string().clone();
            
            for i in (0..(seq.len() - kmer_length as usize + 1)).step_by(1 as usize) {
                let kmer = seq[i..(((i as u32)+kmer_length) as usize)].to_string();
                let kmer_hash = hash(&kmer);
                if hm_index.contains_key(&kmer_hash){
                    sum += hm_index.get(&kmer_hash).unwrap();
                }
            }
        }
        counter = counter+1;
    }
    pb.finish_print(&format!("[{}] alignment completed", Local::now().format("%Y-%m-%d][%H:%M:%S")).green());
    //println!("{} total hashed", sum);

    let elapsed = now.elapsed().as_millis();
    println!("time: {}", elapsed-elapsed_read);
    println!("output: {}", sum);
}

pub fn create_index(input_file: &str, kmer_length: u32) -> IntMap<u32, u64> {
    let mut seq: String;
    let mut counter = 0;
    let mut sum = 0;
    let line_count = 2600000;
    //let hm: IntMap<u64, u64, nohash_hasher::BuildNoHashHasher<u64>> = IntMap::new();
    
    let mut hm_index: IntMap<u32, u64> = IntMap::default();
    
    let mut pb = ProgressBar::new(line_count as u64);
    pb.set_max_refresh_rate(Some(Duration::from_millis(100)));
    pb.message("     -> build index | ");

    let now = Instant::now();
    let f = BufReader::new(File::open(input_file).unwrap());
    for (_i, line) in f.lines().enumerate() {
        
        if counter % 4 == 1 {
            pb.inc();
            seq = line.unwrap().to_string().clone();
            
            for i in (0..(seq.len() - kmer_length as usize + 1)).step_by(1 as usize) {
                let kmer = seq[i..(((i as u32)+kmer_length) as usize)].to_string();
                let kmer_hash = hash(&kmer);
                hm_index.insert(kmer_hash, sum as u64);
                sum += 1;
            }
        }
        counter = counter+1;
        
    }
    pb.finish_print(&format!("[{}] alignment completed", Local::now().format("%Y-%m-%d][%H:%M:%S")).green());
    println!("{} total hashed", sum);
    println!("index keys: {}", hm_index.len());
    let elapsed = now.elapsed().as_millis();

    return hm_index;
}


pub fn hash_kmers(seq: String, hm_index: IntMap<u32, u64>, kmer_length: u32, step_size: u32)-> u64 {
    
    let mut kmer;
    let mut sum = 0;
    for i in (0..(seq.len() - kmer_length as usize + 1)).step_by(step_size as usize) {
        kmer = seq[i..(((i as u32)+kmer_length) as usize)].to_string();
        let kmer_hash = hash(&kmer);
        if hm_index.contains_key(&kmer_hash){
            sum += hm_index.get(&kmer_hash).unwrap();
        }
    }

    return sum as u64
}

pub fn hash(input: &str) -> u32 {
    return city_hash_64(input.as_bytes()) as u32;
}


fn main() {
    //let mut hashtable = IntMap::new();
    let matches = cliparameters::cli();
    let input_file = matches.value_of("input-file").unwrap();
    let kmer_length = matches.value_of("k-mer").unwrap_or("28").parse::<u32>().unwrap();

    let now = Instant::now();
    let hm_index = create_index(input_file, kmer_length);
    let elapsed_index = now.elapsed().as_millis();
    println!("Elapsed time (index): {}", elapsed_index);

    for n in 1..3 {
        let x1 = read_fastq(input_file, hm_index.clone(), kmer_length);
        //println!("{}", format!("[{}] Done! Elapsed time: {}ms", Local::now().format("%Y-%m-%d][%H:%M:%S"), now.elapsed().as_millis()).green());
    }
}
