
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

pub fn read_fastq(input_file: &str, kmer_length: u32, eq_classes: &HashMap<u32, u32>, eq_elements: &HashMap<u32, Vec<u32>>, _kmer_offset: &u32, step_size: u32, sensitivity: u32) -> (i64, HashMap<u32, u32>, i64) {

    let mut transcript_counts: HashMap<u32, u32> = HashMap::new();

    let mut seq: String;
    let mut counter = 0;
    let mut read_length = 0;
    //let mut ncounter = 0;
    
    let mut line_count: i64 = 0;
    let f = BufReader::new(File::open(input_file).unwrap());
    for (_i, line) in f.lines() {
        line_count = line_count+1;

        if line_count == 2 {
            read_length = line.unwrap().to_string().len();
        }
    }

    let mut total_match_count = 0;
    let mut unique_counter = 0;

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
                //ncounter += 1;
                let na_replaced = remove_n(seq.clone());
                for seq in na_replaced {
                    let (match_counter, matches, is_unique) = get_matches(seq, eq_classes, eq_elements, kmer_length, step_size);
                    total_match_count = total_match_count + match_counter;
                    
                    if matches.len() > 0 {
                        if matches.len() == 1 {
                            let temp = matches[0];
                            let counter_tmp = transcript_counts.entry(temp).or_insert(0);
                            *counter_tmp += 1;
                            break;
                        }
                        else {
                            let (top_transcript, max_count) = most_frequent(&matches,1);
                            let temp = top_transcript.clone() as u32;
                            if max_count > 0 {
                                let counter_tmp = transcript_counts.entry(temp).or_insert(0);
                                *counter_tmp += 1;
                            }
                            break;
                        }
                    }
                }
            }
            else {
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
                
                total_match_count = total_match_count + match_counter;
                
                if matches.len() > 0 {
                    if matches.len() == 1 && false {
                        let temp = matches[0];
                        let counter_tmp = transcript_counts.entry(temp).or_insert(0);
                        *counter_tmp += 1;
                    }
                    else {
                        let (top_transcript, max_count) = most_frequent(&matches,1);
                        let temp = top_transcript.clone() as u32;
                        if max_count > sensitivity {
                            let counter_tmp = transcript_counts.entry(temp).or_insert(0);
                            *counter_tmp += 1;
                            if is_unique {
                                unique_counter += 1;
                            }
                        }
                    }
                }
            }
        }

        counter = counter+1;
    }
    pb.finish_print(&format!("[{}] alignment completed", Local::now().format("%Y-%m-%d][%H:%M:%S")).green());
    println!("");
    
    //println!("N counts: {}", ncounter);
    return (line_count, transcript_counts, unique_counter);
}

pub fn get_matches(seq: String, eq_classes: &HashMap<u32, u32>, eq_elements: &HashMap<u32, Vec<u32>>, kmer_length: u32, step_size: u32) -> (u32, Vec<u32>, bool){
    
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
                    let mut u: Vec<u32> = vec![];
                    u.push(*p);
                    matches.extend(u);
                    break;
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
                    let mut u: Vec<u32> = vec![];
                    u.push(*p);
                    matches.extend(u);
                    break;
                }
            }
        }
    }

    return (match_counter, matches, is_unique);
}

pub fn most_frequent<T: std::hash::Hash + std::cmp::Eq + std::cmp::Ord>(array: &[T], n: u32) -> (&T, u32) {

    let mut map = HashMap::new();

    for value in array {
        let counter = map.entry(value).or_insert(0);
        *counter += 1;
    }

    let mut heap: BinaryHeap<_> = map.values().collect();
    let mut max = heap.pop().unwrap();

    for _i in 1..n {
        max = heap.pop().unwrap();
    }

    let top_element = map.iter()
        .find_map(|(key, &val)| if val == *max { Some(key) } else { None })
        .unwrap();
    
    return (top_element, *max);
}

pub fn hash(input: &str) -> u32 {
    let mut hasher = DefaultHasher::new();
    input.hash(&mut hasher);
    return hasher.finish() as u32;
}
