use std::cmp;
use std::fs::File;
use std::iter::FromIterator;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::error::Error;
use rand::Rng;
use pbr::ProgressBar;
use std::time::{Instant};

use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};

use std::collections::HashMap;
use std::collections::HashSet;

use std::collections::{BinaryHeap};

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

pub fn read_fastq(input_file: &str, contig_length: i32, equivalence_classes: &HashMap<u32, Vec<i32>>, contig_offset: &u32) -> (i64, HashMap<u32, u32>) {

    let mut transcript_counts: HashMap<u32, u32> = HashMap::new();

    let mut seq = "".to_owned();
    let mut counter = 0;
    
    let mut line_count: i64 = 0;
    let f = BufReader::new(File::open(input_file).unwrap());
    for (i, line) in f.lines().enumerate() {
        line_count = line_count+1;
    }

    println!("Total sequence count: {}", line_count/4);

    //let mut pb = ProgressBar::new((line_count/4) as u64);
    
    
    //let step_size = contig_length as u32;
    let step_size = 1;
    let mut total_match_count = 0;

    let now = Instant::now();

    let f = BufReader::new(File::open(input_file).unwrap());
    for (i, line) in f.lines().enumerate() {
        
        if counter % 4 == 1 {
            //pb.inc();
            seq = line.unwrap().to_string().clone();
            
            let (match_counter, matches) = get_matches(seq, equivalence_classes, contig_length, step_size);
            total_match_count = total_match_count + match_counter;
            
            if matches.len() > 0 {
                let (top_transcript, max_count) = most_frequent(&matches,1);
                let temp = top_transcript.clone() as u32;
                if max_count > 3 {
                    let counter_tmp = transcript_counts.entry(temp).or_insert(0);
                    *counter_tmp += 1;
                }
            }
        }

        counter = counter+1;
    }

    println!("\nMatched sequences: {}", total_match_count);

    return (line_count, transcript_counts);
}

pub fn get_matches(seq: String, equivalence_classes: &HashMap<u32, Vec<i32>>, contig_length: i32, step_size: u32) -> (i32, Vec<i32>){
    
    let mut matches: Vec<i32> = vec![];
    let mut hasher = DefaultHasher::new();
    let mut rev_contig = "".to_string();
    let mut contig = "".to_string();
    let mut match_counter = 0;
    let mut u = &vec![];
    let mut hash_result = 0;

    for i in (0..((seq.len() as i32 - contig_length + 1) as usize)).step_by(step_size as usize) {
        
        contig = seq[i..(((i as i32)+contig_length) as usize)].to_string();
        
        hasher = DefaultHasher::new();
        contig.hash(&mut hasher);
        hash_result = hasher.finish() as u32;
        if equivalence_classes.contains_key(&hash_result) {
            match_counter = match_counter + 1;
            u = equivalence_classes.get(&hash_result).unwrap();
            matches.extend(u);
        }
        else{
            rev_contig = reverse_read(contig);
            hasher = DefaultHasher::new();
            rev_contig.hash(&mut hasher);
            hash_result = hasher.finish() as u32;
            if equivalence_classes.contains_key(&hash_result) {
                match_counter = match_counter + 1;
                u = equivalence_classes.get(&hash_result).unwrap();
                matches.extend(u);
            }
        }
    }

    return (match_counter, matches);
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