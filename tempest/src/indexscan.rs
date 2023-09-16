use std::collections::HashMap;
use std::fs::File;
use std::io::{Error, Write};


use nohash_hasher::IntMap;

pub fn index_stats(transcripts: &Vec<String>, eq_classes: &IntMap<u32, u32>, eq_elements: &IntMap<u32, Vec<u32>>, transcript_kmers: &IntMap<u32, u32>, kmer_length: &u32)  -> std::io::Result<()>{
    println!("IS transcripts: {} | eq_classes: {} | eq_elements: {} | k-mer length: {}", transcripts.len(), eq_classes.len(), eq_elements.len(), kmer_length);

    let mut unique = 0;
    let mut ambiguous = 0;
    let mut counts = vec![0;300];
    for k in eq_classes.keys() {
        if eq_classes.get(k).unwrap() < &(transcripts.len() as u32) {
            unique += 1;
        }
        else {
            ambiguous += 1;
        }
    }

    let mut once = 0;
    for k in eq_elements.keys() {
        if eq_elements.get(k).unwrap().len() < counts.len() {
            counts[eq_elements.get(k).unwrap().len()] += 1;
        }
    }

    println!("Uniques: {} | Ambiguous: {}", unique, ambiguous);

    let mut unique_count = vec![0;1000];
    let mut temp_map = IntMap::default();
    for k in eq_classes.keys() { 
        if (eq_classes.get(k).unwrap().clone() as usize) < transcripts.len() {
            *temp_map.entry(eq_classes.get(k).unwrap().clone()).or_insert(0) += 1;
        }
    }
    
    let mut output = File::create("output/transcript_info_40.tsv")?;
    for k in transcript_kmers.keys() {
        if temp_map.contains_key(k) {
            write!(output, "{}\t{}\t{}\n", transcripts[k.clone() as usize], temp_map.get(k).unwrap(), transcript_kmers.get(k).unwrap())?;
        }
        else {
            write!(output, "{}\t{}\t{}\n", transcripts[k.clone() as usize], 0, transcript_kmers.get(k).unwrap())?;
        }
    }

    for k in temp_map.keys() {
        let tt = temp_map.get(k).unwrap().clone() as usize;
        if tt < unique_count.len() {
            unique_count[tt] += 1;
        }
    }

    println!("Total identifiable: {}: ", temp_map.len());
    //println!("Unique seg: {:?}", unique_count);

    Ok(())
}