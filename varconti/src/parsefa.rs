
use std::fs::File;
use std::io::{BufRead, BufReader};

use pbr::ProgressBar;
use std::time::Duration;

use std::hash::{Hash, Hasher};
use std::collections::HashMap;
use std::collections::hash_map::DefaultHasher;

extern crate chrono;
use chrono::Local;
use colored::*;

pub fn read_fa(input_file: &str, kmer_length: u32) -> (Vec<String>, HashMap<u32, u32>, HashMap<u32, Vec<u32>>) {

    let mut eq_classes: HashMap<u32, u32> = HashMap::new();
    let mut eq_elements: HashMap<u32, Vec<u32>> = HashMap::new();
    let mut all = 0;
    
    let mut sequence_count = 0;
    let f = BufReader::new(File::open(input_file).unwrap());
    for line in f.lines() {

        let lw = line.unwrap();
        if &lw[..1] == ">" {
            sequence_count += 1;
        }
    }

    let mut pb = ProgressBar::new(sequence_count as u64);
    pb.set_max_refresh_rate(Some(Duration::from_millis(200)));
    pb.message("     -> k-mer processing | ");
    
    let mut transcripts: Vec<String> = vec![];
    let mut seq = "".to_owned();
    let mut first = true;
    let mut total = 0;

    let mut counter: u32 = 0;
    let f = BufReader::new(File::open(input_file).unwrap());
    let mut tc = 0;

    for line in f.lines() {
        tc += 1;
        let lw = line.unwrap();

        if &lw[..1] == ">" {
            
            let v: Vec<&str> = lw.split(" ").collect();
            let ensembl_transcript_id = v[0].to_string().replace(">", "");
            transcripts.push(ensembl_transcript_id);

            if !first {
                pb.set(counter as u64);
                
                let kmers = split_kmers(&seq, kmer_length);

                for j in 0..kmers.len() {
                    total += 1;
                    let kmer_hash = hash(&kmers[j]);
                    if eq_classes.contains_key(&kmer_hash) {
                        let elements_pointer = eq_classes.get(&kmer_hash).unwrap();
                        if eq_elements.contains_key(&elements_pointer) {
                            // This means there is a set larger than 1
                            // 1. Create a copy of the vec and add sequence number
                            // 2. Calculate hash of vec
                            // 3. Update entry for contig hash in eq_classes
                            // 4. Insert new set in eq_classes_map if does not exist
                            let mut elements = eq_elements.get_mut(&elements_pointer).unwrap().clone();
                            
                            if elements.len() < 30 {
                                let mut testin = false;
                                for k in 0..elements.len() {
                                    if elements[k] == counter{
                                        testin = true;
                                        break;
                                    }
                                }
                                if !testin {
                                    elements.push(counter);
                                    let elements_hash = listhash(&elements);
                                    eq_classes.insert(kmer_hash, elements_hash);
                                    if !eq_elements.contains_key(&elements_hash) {
                                        eq_elements.insert(elements_hash, elements);
                                    }
                                }
                            }
                        }
                        else {
                            // new EQ class
                            let mut elements = vec![];
                            elements.push(elements_pointer.clone());
                            elements.push(counter.clone());
                            
                            let elements_hash: u32 = listhash(&elements);
                            eq_classes.insert(kmer_hash, elements_hash);
                            
                            if !eq_elements.contains_key(&elements_hash) {
                                eq_elements.insert(elements_hash, elements.clone());
                            }
                        }
                    }
                    else {
                        eq_classes.insert(kmer_hash, counter);
                    }
                    all = all+1;
                }
                counter += 1;
            }

            first = false;
            seq = "".to_owned();
        }
        else {
            seq = seq.to_owned();
            seq.push_str(&lw);
        }
    }

    pb.finish();

    let mut pb = ProgressBar::new(eq_classes.len() as u64);
    pb.set_max_refresh_rate(Some(Duration::from_millis(200)));
    pb.message("     -> cleaning classes | ");
    //let mut total = 0;
    let mut counter = 1;
    for _k in eq_classes.keys() {
        counter = counter + 1;
        if counter % 100000 == 0 {
            pb.set(counter);
        }
        //if eq_classes_map.contains_key(eq_classes.get(k).unwrap()) {
        //    total += 1;
        //}
    }

    let mut sum_elements = 0;
    for k in eq_elements.keys() {
        sum_elements += eq_elements.get(k).unwrap().len();
    }


    pb.finish_print(&format!("[{}] k-mer hash completed", Local::now().format("%Y-%m-%d][%H:%M:%S")).green());
    println!("");
    println!("Total: {}", total);
    println!("EQ elements: {}", eq_elements.len());
    println!("Avg length elements: {}", sum_elements/eq_elements.len());
    return (transcripts, eq_classes, eq_elements);
}

pub fn split_kmers(sequence: &str, kmer_length: u32) -> Vec<String> {
    let mut kmers = vec![];
    if sequence.len() >= kmer_length as usize {
        for i in 0..((sequence.len() as u32 - kmer_length + 1) as usize) {
            let kmer = sequence[i..(((i as u32)+kmer_length) as usize)].to_string();
            kmers.push(kmer);
        }
    }
    return kmers;
}

pub fn hash(input: &str) -> u32 {
    let mut hasher = DefaultHasher::new();
    input.hash(&mut hasher);
    return hasher.finish() as u32;
}

pub fn listhash(input: &Vec<u32>) -> u32 {
    let mut h = DefaultHasher::new();
    Hash::hash_slice(input, &mut h);
    return h.finish() as u32;
}