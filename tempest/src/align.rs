
use std::fs::File;
use std::io::{BufRead, BufReader};
use pbr::ProgressBar;

use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};
use rustc_hash::FxHashMap;
use std::collections::HashMap;

use std::collections::{BinaryHeap};

extern crate chrono;
use chrono::Local;
use colored::*;
use std::time::Duration;

use std::sync::{Arc, Mutex, RwLock};

use fxhash::FxHasher32;
use fxhash;

use nohash_hasher::IntMap;
use cityhash::cityhash_1_1_1::city_hash_64;

use rayon::prelude::*;

pub fn reverse_read(read: String) -> String {
    let s:String = read.chars()
        .map(|x| match x { 
            'A' => 'T', 
            'T' => 'A',
            'G' => 'C', 
            'C' => 'G',
            _ => x
    }).rev().collect();
    return s;
}

pub fn remove_n(read: String) -> Vec<String> {
    let t:String = read.chars()
        .map(|x| match x { 'N' => 'T', _ => x}).collect();

    let a:String = read.chars()
        .map(|x| match x { 'N' => 'A', _ => x}).collect();

    let c:String = read.chars()
        .map(|x| match x { 'N' => 'C', _ => x}).collect();

    let g:String = read.chars()
        .map(|x| match x { 'N' => 'G', _ => x}).collect();

    let mut result: Vec<String> = vec![];
    result.push(t);
    result.push(a);
    result.push(g);
    result.push(c);

    return result;
}

pub fn read_fastq(input_file: &str, kmer_length: usize, eq_classes: &FxHashMap<u64, u32>, eq_elements: &FxHashMap<u32, Vec<u32>>, _kmer_offset: &u32, step_size: u32, sensitivity: u32) -> (i64, Vec<Vec<u32>>) {

    let mut transcript_counts_unique: FxHashMap<u32, u32> = FxHashMap::default();
    let mut transcript_counts: FxHashMap<u32, u32> = FxHashMap::default();
    let mut alignment_matches: Vec<Vec<u32>> = vec![];
    
    let mut seq: String;
    let mut counter = 0;
    let mut read_length = 0;
    //let mut ncounter = 0;
    
    let mut line_count: i64 = 0;
    let f = BufReader::new(File::open(input_file).unwrap());
    for line in f.lines() {
        line_count += 1;

        if line_count == 2 {
            read_length = line.unwrap().to_string().len();
        }
    }

    let mut total_match_count = 0;
    let mut unique_counter = 0;

    let mut max_tc = 0;

    //let mut pb = ProgressBar::new(line_count as u64/4);
    //pb.set_max_refresh_rate(Some(Duration::from_millis(200)));
    //pb.message("     -> aligning | ");

    let pb = Arc::new(Mutex::new(ProgressBar::new((line_count / 4) as u64)));
    {
        let pb = Arc::clone(&pb);
        pb.lock().unwrap().set_max_refresh_rate(Some(Duration::from_millis(200)));
        pb.lock().unwrap().message("    -> aligning | ");
    }

    let mut align_strategy = vec![0; read_length];

    let f = BufReader::new(File::open(input_file).unwrap());
    let lines: Vec<_> = f.lines().collect();
    
    let alignment_matches: Vec<Vec<u32>> = lines.into_par_iter()
        .enumerate()
        .filter_map(|(counter, line)| {
            if counter % 4 == 1 {
                let seq = line.unwrap().to_string();

                let mut alignments: Vec<Vec<u32>> = Vec::new();

                if seq.contains("N") {
                    let na_replaced = remove_n(seq);
                    for seq in na_replaced {
                        process_seq(seq, eq_classes, eq_elements, kmer_length, sensitivity, &mut alignments);
                    }
                } else {
                    process_seq(seq, eq_classes, eq_elements, kmer_length, sensitivity, &mut alignments);
                }
                pb.lock().unwrap().inc();
                return Some(alignments);
            }
            None
        })
        .flatten()
        .collect();

    pb.lock().unwrap().finish_print(&format!("[{}] alignment completed", Local::now().format("%Y-%m-%d][%H:%M:%S")).green());
    println!("");

    return (line_count, alignment_matches)
}


fn process_seq(seq: String, eq_classes: &FxHashMap<u64, u32>, eq_elements: &FxHashMap<u32, Vec<u32>>, kmer_length: usize, sensitivity: u32, alignment_matches: &mut Vec<Vec<u32>>) {
    let (match_counter, matches, is_unique) = get_matches_control(seq, eq_classes, eq_elements, kmer_length, kmer_length, sensitivity);

    if matches.len() > 0 {
        if matches.len() == 1 {
            alignment_matches.push(matches);
        } else {
            alignment_matches.push(most_frequent(&matches));
        }
    }
}

pub fn get_matches_control(seq: String, eq_classes: &FxHashMap<u64, u32>, eq_elements: &FxHashMap<u32, Vec<u32>>, kmer_length: usize, step_size: usize, sensitivity: u32) -> (u32, Vec<u32>, bool) {
    let (mut match_counter, mut matches, mut is_unique) = get_matches(seq.clone(), eq_classes, eq_elements, kmer_length, 1);
    if match_counter < sensitivity {
        let (match_counter_t, matches_t, is_unique_t) = get_matches(seq.clone(), eq_classes, eq_elements, kmer_length, kmer_length/2);
        match_counter = match_counter_t;
        matches = matches_t;
        is_unique = is_unique_t;
    }
    if match_counter < sensitivity {
        let (match_counter_t, matches_t, is_unique_t) = get_matches(seq.clone(), eq_classes, eq_elements, kmer_length, kmer_length/4);
        match_counter = match_counter_t;
        matches = matches_t;
        is_unique = is_unique_t;
    }
    if match_counter < sensitivity {
        let (match_counter_t, matches_t, is_unique_t) = get_matches(seq, eq_classes, eq_elements, kmer_length, step_size);
        match_counter = match_counter_t;
        matches = matches_t;
        is_unique = is_unique_t;
    }

    return (match_counter, matches, is_unique);
}

pub fn get_matches(seq: String, eq_classes: &FxHashMap<u64, u32>, eq_elements: &FxHashMap<u32, Vec<u32>>, kmer_length: usize, step_size: usize) -> (u32, Vec<u32>, bool){
    
    let mut matches: Vec<u32> = vec![];
    let mut rev_kmer;
    let mut kmer;
    let mut match_counter = 0;

    let mut rev_seq = String::new();
    let mut reverse_mode = false;
    let mut pos_mode = false;

    let mut is_unique = true;

    for i in (0..(seq.len() - kmer_length + 1)).step_by(step_size) {
        
        // Positive mode
        if !reverse_mode {
            kmer = &seq[i..i + kmer_length];
            
            let kmer_hash = hash(kmer);
            if let Some(p) = eq_classes.get(&kmer_hash) {
                pos_mode = true;
                match_counter += 1;
                if let Some(elements) = eq_elements.get(p) {
                    is_unique = false;
                    matches.extend(elements);
                } else {
                    matches = vec![*p];
                    is_unique = true;
                    break;
                }
            }
        }
        
        // Reverse mode
        if !pos_mode {
            if i == 0 {
                rev_seq = reverse_read(seq.clone());
            }
            let index = rev_seq.len() - kmer_length - i;
            rev_kmer = &rev_seq[index..index + kmer_length];
            
            let kmer_hash = hash(rev_kmer);
            if let Some(p) = eq_classes.get(&kmer_hash) {
                reverse_mode = true;
                match_counter += 1;
                if let Some(elements) = eq_elements.get(p) {
                    is_unique = false;
                    matches.extend(elements);
                } else {
                    matches = vec![*p];
                    is_unique = true;
                    break;
                }
            }
        }
    }

    (match_counter, matches, is_unique)
}

fn most_frequent(array: &Vec<u32>) -> Vec<u32> {
    let mut map:FxHashMap<u32, u32> = FxHashMap::default();
    for value in array {
        *map.entry(*value).or_insert(0) += 1;
    }

    let max_count = map.values().cloned().max().unwrap_or(0);
    map.into_iter()
        .filter(|&(_, count)| count == max_count)
        .map(|(num, _)| num)
        .collect()
}

fn most_frequent2(array: &Vec<u32>) -> Vec<u32> {
    let mut map:FxHashMap<u32, u32> = FxHashMap::default();
    for value in array {
        *map.entry(*value).or_insert(0) += 1;
    }
    let heap: BinaryHeap<_> = map.values().collect();
    let max = heap.peek().unwrap();
    let top_elements = find_keys_for_value(&map, max).to_vec();
    let mut result = vec![];
    for e in top_elements {
        result.push(e.clone());
    }
    //result.sort();
    return result;
}

fn find_keys_for_value<'a>(map: &'a FxHashMap<u32, u32>, value: &'a u32) -> Vec<&'a u32> {
    map.iter()
        .filter_map(|(key, &val)| if &val == value { Some(key) } else { None })
        .collect()
}

pub fn hash(input: &str) -> u64 {
    return city_hash_64(input.as_bytes());
}
