
use std::fs::File;
use std::io::{BufRead, BufReader};

use pbr::ProgressBar;
use std::time::Duration;

use std::hash::{Hash, Hasher};
use std::collections::{HashMap, HashSet};
use std::collections::hash_map::DefaultHasher;
use rustc_hash::FxHashMap;

extern crate chrono;
use chrono::Local;
use colored::*;

use fxhash::FxHasher;
use fxhash;

use nohash_hasher::IntMap;
use cityhash::cityhash_1_1_1::city_hash_64;

pub fn read_fa(input_file: &str, kmer_length: usize) -> (Vec<String>, FxHashMap<u64, u32>, FxHashMap<u32, Vec<u32>>, FxHashMap<u32, u32>, Vec<usize>) {
    let mut transcript_length: Vec<usize> = vec![];
    
    let mut transcript_kmers: FxHashMap<u32, u32> = FxHashMap::default();
    let mut eq_classes: FxHashMap<u64, u32> = FxHashMap::default();
    let mut eq_elements: FxHashMap<u32, Vec<u32>> = FxHashMap::default();

    let sequence_count = BufReader::new(File::open(input_file).unwrap())
        .lines()
        .filter_map(Result::ok)
        .filter(|lw| lw.starts_with('>'))
        .count();

    let mut pb = ProgressBar::new(sequence_count as u64);
    pb.set_max_refresh_rate(Some(Duration::from_millis(200)));
    pb.message("     -> k-mer processing | ");
    
    let mut transcripts: Vec<String> = vec![];
    let mut seq = String::new();
    let mut counter: u32 = 0;

    let f = BufReader::new(File::open(input_file).unwrap());
    
    for line in f.lines() {
        let lw = line.unwrap().to_uppercase();

        if lw.starts_with('>') {
            let ensembl_transcript_id = lw.split_whitespace().next().unwrap().replace(">", "");
            transcripts.push(ensembl_transcript_id);

            if !seq.is_empty() {
                pb.set(counter as u64);
                process_sequence(&mut seq, kmer_length, counter, &mut transcript_length, &mut eq_classes, &mut eq_elements, &mut transcript_kmers);
                counter += 1;
            }
            seq.clear();
        }
        else {
            seq.push_str(&lw);
        }
    }

    if !seq.is_empty() {
        pb.set(counter as u64);
        process_sequence(&mut seq, kmer_length, counter, &mut transcript_length, &mut eq_classes, &mut eq_elements, &mut transcript_kmers);
    }

    pb.finish();

    let mut pb = ProgressBar::new(eq_classes.len() as u64);
    pb.set_max_refresh_rate(Some(Duration::from_millis(200)));
    pb.message("     -> cleaning classes | ");
    
    let (eq_classes, eq_elements) = purge_disconnected(eq_classes, eq_elements);

    let sum_elements: usize = eq_elements.values().map(|v| v.len()).sum();
    
    pb.finish_print(&format!("[{}] k-mer hash completed", Local::now().format("%Y-%m-%d][%H:%M:%S")).green());
    println!("");
    println!("Total: {}", transcript_kmers.values().map(|v| *v as usize).sum::<usize>());
    println!("EQ classes: {}", eq_classes.len());
    println!("EQ elements: {}", eq_elements.len());
    println!("Avg length elements: {}", sum_elements/eq_elements.len());
    println!("final counter (max element): {}", counter);
    return (transcripts, eq_classes, eq_elements, transcript_kmers, transcript_length);
}

fn process_sequence(seq: &mut String, kmer_length: usize, counter: u32, transcript_length: &mut Vec<usize>, eq_classes: &mut FxHashMap<u64, u32>, eq_elements: &mut FxHashMap<u32, Vec<u32>>, transcript_kmers: &mut FxHashMap<u32, u32>) {
    let kmers = split_kmers(&seq, kmer_length);
    transcript_length.push(seq.len());

    for j in 0..kmers.len() {
        let kmer_hash = hash(&kmers[j]);
        if let Some(elements_pointer) = eq_classes.get(&kmer_hash) {
            if let Some(elements) = eq_elements.get_mut(&elements_pointer) {
                // This means there is a set larger than 1
                // 1. Create a copy of the vec and add sequence number
                // 2. Calculate hash of vec
                // 3. Update entry for contig hash in eq_classes
                // 4. Insert new set in eq_classes_map if does not exist
                let mut elements = elements.clone();
                
                // test if element already in eqc
                let testin = elements.contains(&counter);
                if !testin {
                    elements.push(counter);
                    elements.sort_unstable();
                    let elements_hash = listhash(elements.clone());
                    eq_classes.insert(kmer_hash, elements_hash);
                    if !eq_elements.contains_key(&elements_hash) {
                        eq_elements.insert(elements_hash, elements);
                    }
                }
            }
            else {
                // new EQ class
                let mut elements = vec![];
                elements.push(elements_pointer.clone());
                elements.push(counter.clone());
                elements.sort_unstable();

                let elements_hash: u32 = listhash(elements.clone());
                eq_classes.insert(kmer_hash, elements_hash);
                
                if !eq_elements.contains_key(&elements_hash) {
                    eq_elements.insert(elements_hash, elements.clone());
                }
            }
        }
        else {
            eq_classes.insert(kmer_hash, counter.clone());
        }
    }

    transcript_kmers.insert(counter, kmers.len() as u32);
}

pub fn split_kmers(sequence: &str, kmer_length: usize) -> Vec<String> {
    let mut kmers = Vec::new();
    let seq_len = sequence.len();
    
    if seq_len >= kmer_length {
        for i in 0..=(seq_len - kmer_length) {
            let kmer = sequence[i..(i + kmer_length)].to_string();
            kmers.push(kmer);
        }
    }   
    kmers
}

pub fn split_kmers_old(sequence: &str, kmer_length: u32) -> Vec<String> {
    let mut kmers = vec![];
    if sequence.len() >= kmer_length as usize {
        for i in 0..((sequence.len() as u32 - kmer_length + 1) as usize) {
            let kmer = sequence[i..(((i as u32)+kmer_length) as usize)].to_string();
            kmers.push(kmer);
        }
    }
    return kmers;
}

pub fn purge_disconnected(eq_classes: FxHashMap<u64, u32>, eq_elements: FxHashMap<u32, Vec<u32>>) -> (FxHashMap<u64, u32>, FxHashMap<u32, Vec<u32>>){
    
    let mut hh = HashSet::new();
    for k in eq_classes.keys() {
        hh.insert(eq_classes.get(k).unwrap());
    }
    
    let mut eqtemp = eq_elements.clone();
    for k in eq_elements.keys() {
        if !hh.contains(&k) {
            eqtemp.remove(k);
        }
    }
    return (eq_classes, eqtemp);
}


pub fn hash(input: &str) -> u64 {
    return city_hash_64(input.as_bytes());
}

pub fn listhash_old(input: Vec<u32>) -> u32 {
    let s = format!("{:?}", input);
    let mut hasher = FxHasher::default();
    hasher.write_u32(0);
    &s.hash(&mut hasher);
    return hasher.finish() as u32;
}

pub fn listhash(input: Vec<u32>) -> u32 {
    let mut hasher = FxHasher::default();
    hasher.write_u32(0);
    for &number in &input {
        hasher.write_u32(number);
    }
    hasher.finish() as u32
}