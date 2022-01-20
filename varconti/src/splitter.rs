use std::cmp;
use std::io::{BufRead, BufReader, BufWriter, Write};

use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};

use std::collections::HashMap;
use std::collections::HashSet;

use serde_json::{Result, Value};

use std::fs::OpenOptions;
use std::fs;
use std::error::Error;

use std::fs::File;
use serde_json;

extern crate flate2;

use flate2::write::ZlibEncoder;
use flate2::read::ZlibDecoder;
use flate2::Compression;

extern crate rustc_serialize;
extern crate bincode;

use std::sync::{Arc, RwLock};
use threadpool::ThreadPool;
use std::sync::mpsc::channel;
use pbr::ProgressBar;

use fnv::FnvHashMap;

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

pub fn split_contigs(sequence: &str, contig_length: i32) -> Vec<String> {

    let mut contigs = vec![];
    if sequence.len() >= contig_length as usize {
        for i in 0..((sequence.len() as i32 - contig_length + 1) as usize) {
            let contig = sequence[i..(((i as i32)+contig_length) as usize)].to_string();
            contigs.push(contig);
        }
    }

    return contigs;
}

pub fn contig_index(sequences: Vec<String>, contig_length: i32)  -> std::io::Result<()> {
    
    let mut hash_result = 0;
    //let mut equivalence_classes: HashMap<u32, HashSet<i32>> = HashMap::new();
    let mut equivalence_classes: FnvHashMap<u32, HashSet<i32>> = FnvHashMap::default();
    let mut uni = 0;
    let mut all = 0;
    
    let mut pb = ProgressBar::new(sequences.len() as u64);

    for i in 1..sequences.len() {
        pb.inc();
        let contigs = split_contigs(&sequences[i], contig_length);
        for j in 1..contigs.len() {
            let mut hasher = DefaultHasher::new();
            contigs[j].hash(&mut hasher);
            hash_result = hasher.finish() as u32;
            if equivalence_classes.contains_key(&hash_result) {
                let mut temp = equivalence_classes.get_mut(&hash_result).unwrap();
                temp.insert(i as i32);
            }
            else {
                uni = uni + 1;
                let mut temp = HashSet::new();
                temp.insert(i as i32);
                equivalence_classes.insert(hash_result, temp);
            }
            all = all+1;
        }

        let contigs = split_contigs(&sequences[i], 12);
        for j in 1..contigs.len() {
            let mut hasher = DefaultHasher::new();
            contigs[j].hash(&mut hasher);
            hash_result = hasher.finish() as u32;
            if equivalence_classes.contains_key(&hash_result) {
                let mut temp = equivalence_classes.get_mut(&hash_result).unwrap();
                temp.insert(i as i32);
            }
            else {
                uni = uni + 1;
                let mut temp = HashSet::new();
                temp.insert((i+1000000) as i32);
                equivalence_classes.insert(hash_result, temp);
            }
            all = all+1;
        }
    }
    pb.finish_print("Completed index generation.");

    println!("unique: {:.2}", uni as f64 / all as f64);
    //let res = bincode::serialize(&equivalence_classes).unwrap();
    //fs::write("binary.idx", res).unwrap();

    //let file = File::create("index.json")?;
    //let mut writer = BufWriter::new(file);
    //serde_json::to_writer(&mut writer, &equivalence_classes)?;
    //writer.flush()?;
    Ok(())
}


pub fn contig_index_vec(sequences: Vec<String>, contig_length: i32) {
    
    let mut hash_result = 0;
    let mut equivalence_classes: HashMap<u32, Vec<i32>> = HashMap::new();
    //let mut equivalence_classes: HashMap<u32, Vec<i32>> = HashMap::with_capacity(4000000);
    //let mut equivalence_classes: FnvHashMap<u32, Vec<i32>> = FnvHashMap::default();
    let mut uni = 0;
    let mut all = 0;
    
    let mut pb = ProgressBar::new(sequences.len() as u64);

    for i in 0..sequences.len() {
        pb.inc();
        let contigs = split_contigs(&sequences[i], contig_length);
        for j in 0..contigs.len() {
            let mut hasher = DefaultHasher::new();
            contigs[j].hash(&mut hasher);
            hash_result = hasher.finish() as u32;
            if equivalence_classes.contains_key(&hash_result) {
                let mut temp = equivalence_classes.get_mut(&hash_result).unwrap();
                temp.push(i as i32);
            }
            else {
                uni = uni + 1;
                let mut temp = vec![];
                temp.push(i as i32);
                equivalence_classes.insert(hash_result, temp);
            }
            all = all+1;
        }
    }

    pb.finish_print("Completed index generation.");
    println!("");
    println!("contigs: {} | unique contigs: {} |  % unique: {:.2}", all, uni, uni as f64 / all as f64);
    let res = bincode::serialize(&equivalence_classes).unwrap();
    fs::write("binary.idx", res).unwrap();

    std::thread::spawn(move || drop(sequences));
}
