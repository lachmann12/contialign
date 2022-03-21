
use std::fs::File;
use std::io::{BufRead, BufReader};
use pbr::ProgressBar;

use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};

use std::collections::HashMap;

use std::collections::{BinaryHeap};

extern crate chrono;
use chrono::Local;
use colored::*;
use std::time::Duration;

use fxhash::FxHasher32;
use fxhash;

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

pub fn read_fastq(input_file: &str, kmer_length: u32, eq_classes: &HashMap<u64, u32>, eq_elements: &HashMap<u32, Vec<u32>>, _kmer_offset: &u32, step_size: u32, sensitivity: u32) -> (i64, Vec<Vec<u32>>) {

    let mut transcript_counts_unique: HashMap<u32, u32> = HashMap::new();
    let mut transcript_counts: HashMap<u32, u32> = HashMap::new();
    let mut alignment_matches: Vec<Vec<u32>> = vec![];

    let mut seq: String;
    let mut counter = 0;
    let mut read_length = 0;
    //let mut ncounter = 0;
    
    let mut line_count: i64 = 0;
    let f = BufReader::new(File::open(input_file).unwrap());
    for line in f.lines() {
        line_count = line_count+1;

        if line_count == 2 {
            read_length = line.unwrap().to_string().len();
        }
    }

    let mut total_match_count = 0;
    let mut unique_counter = 0;

    let mut max_tc = 0;

    let mut pb = ProgressBar::new(line_count as u64/4);
    pb.set_max_refresh_rate(Some(Duration::from_millis(200)));
    pb.message("     -> aligning | ");

    let mut align_strategy = vec![0; read_length];

    let f = BufReader::new(File::open(input_file).unwrap());
    for (_i, line) in f.lines().enumerate() {
        
        if counter % 4 == 1 {
            pb.inc();
            seq = line.unwrap().to_string().clone();
            
            if seq.contains("N") {
                
                let na_replaced = remove_n(seq.clone());
                for seq in na_replaced {
                    
                    let (mut match_counter, mut matches, mut is_unique) = get_matches_control(seq.clone(), eq_classes, eq_elements, kmer_length, kmer_length, sensitivity);

                    if matches.len() > 0 {
                        if matches.len() == 1 {
                            alignment_matches.push(matches);
                        }
                        else {
                            let (top_transcripts, max_count) = most_frequent(&matches);
                            alignment_matches.push(top_transcripts);
                        }
                        break;
                    }
                }
            }
            else {
                let (mut match_counter, mut matches, mut is_unique) = get_matches_control(seq.clone(), eq_classes, eq_elements, kmer_length, kmer_length, sensitivity);
                
                if matches.len() > 0 {
                    if matches.len() == 1 {
                        alignment_matches.push(matches);
                    }
                    else {
                        let (top_transcripts, max_count) = most_frequent(&matches);
                        alignment_matches.push(top_transcripts);
                    }
                }
            }
        }
        counter = counter+1;
    }
    pb.finish_print(&format!("[{}] alignment completed", Local::now().format("%Y-%m-%d][%H:%M:%S")).green());
    println!("");

    //return (line_count, transcript_counts_unique, transcript_counts, unique_counter);
    return (line_count, alignment_matches)
}

pub fn get_matches_control(seq: String, eq_classes: &HashMap<u64, u32>, eq_elements: &HashMap<u32, Vec<u32>>, kmer_length: u32, step_size: u32, sensitivity: u32) -> (u32, Vec<u32>, bool) {
    let (mut match_counter, mut matches, mut is_unique) = get_matches(seq.clone(), eq_classes, eq_elements, kmer_length, kmer_length);
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

pub fn get_matches(seq: String, eq_classes: &HashMap<u64, u32>, eq_elements: &HashMap<u32, Vec<u32>>, kmer_length: u32, step_size: u32) -> (u32, Vec<u32>, bool){
    
    let mut matches: Vec<u32> = vec![];
    let mut rev_kmer;
    let mut kmer;
    let mut match_counter = 0;

    let mut rev_seq = "".to_string();
    let mut reverse_mode = false;
    let mut pos_mode = false;

    let mut is_unique = true;

    for i in (0..(seq.len() - kmer_length as usize + 1)).step_by(step_size as usize) {
        
        if !reverse_mode {
            kmer = seq[i..(((i as u32)+kmer_length) as usize)].to_string();
            
            let kmer_hash = hash(&kmer);
            if eq_classes.contains_key(&kmer_hash) {
                pos_mode = true;
                match_counter = match_counter + 1;
                let p = eq_classes.get(&kmer_hash).unwrap();
                if eq_elements.contains_key(p) {
                    is_unique = false;
                    matches.extend(eq_elements.get(p).unwrap());
                }
                else{
                    matches = vec![];
                    matches.push(*p);
                    is_unique = true;
                    //break;
                }
            }
        }
        
        if !pos_mode {
            if i == 0 {
                rev_seq = reverse_read(seq.clone());
            }
            rev_kmer = rev_seq[(rev_seq.len() - kmer_length as usize - i)..(rev_seq.len()-i)].to_string();
            
            let kmer_hash = hash(&rev_kmer);
            if eq_classes.contains_key(&kmer_hash) {
                reverse_mode = true;
                match_counter = match_counter + 1;
                let p = eq_classes.get(&kmer_hash).unwrap();
                if eq_elements.contains_key(p) {
                    is_unique = false;
                    matches.extend(eq_elements.get(p).unwrap());
                }
                else{
                    matches = vec![];
                    matches.push(*p);
                    is_unique = true;
                    //break;
                }
            }
        }
    }

    return (match_counter, matches, is_unique);
}

fn most_frequent(array: &Vec<u32>) -> (Vec<u32>, u32) {
    let mut map:HashMap<u32, u32> = HashMap::new();
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
    return (result, **max);
}

fn find_keys_for_value<'a>(map: &'a HashMap<u32, u32>, value: &'a u32) -> Vec<&'a u32> {
    map.iter()
        .filter_map(|(key, &val)| if &val == value { Some(key) } else { None })
        .collect()
}

pub fn oldhash(input: &str) -> u64 {
    let mut hasher = DefaultHasher::default();
    //hasher.write_u32(0);
    input.hash(&mut hasher);
    return hasher.finish();
    //return fxhash::hash64(input);
}


use std::cmp::PartialEq;
use std::ops::BitXor;
use std::ops::Shl;

pub fn hash(input: &str) -> u64 {
    let mut bits: Vec<u8> = vec![0; 64];
    let mut counter: usize = 0;
    for c in input.chars() {
        //A = 00
        //T = 01
        //C = 10
        //G = 11
        if c == 'T' {
            bits[counter+1] = 1; 
        }
        else if c == 'C' {
            bits[counter] = 1;
        }
        else if c == 'G' {
            bits[counter] = 1;
            bits[counter+1] = 1;
        }
        counter += 2;
    }

    return convert(&bits).unwrap();
}

#[derive(Debug)]
pub enum ConversionError {
    Overflow,
    NonBinaryInput,
}

pub fn convert<T: PartialEq + From<u8> + BitXor<Output=T> + Shl<Output=T> + Clone>(
    bits: &[u8],
) -> Result<T, ConversionError> {
    if bits.len() > (std::mem::size_of::<T>() * 8) {
        return Err(ConversionError::Overflow);
    }
    if bits.iter()
        .filter(|&&bit| bit != 0 && bit != 1).count() > 0 {
        return Err(ConversionError::NonBinaryInput);
    }

    Ok(bits.iter()
        .fold(T::from(0), |result, &bit| {
            (result << T::from(1)) ^ T::from(bit)
        }))
}