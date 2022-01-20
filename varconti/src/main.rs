use std::time::{Instant};

use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};

use std::collections::HashMap;
use std::collections::HashSet;

use std::{thread, time};

mod parsefa;
mod splitter;
mod align;

use std::fs::File;
use std::io::Read;
use std::io::{Error, Write};

extern crate bincode;
use std::iter::FromIterator;
use rand::Rng;



// cargo build --release
// "target/release/varconti.exe"

fn main() -> Result<(), Error> {

    println!("Start");
    
    let mut hasher = DefaultHasher::new();
    let mut contig_length = 22;
    let mut contig_offset = 20;

    let (transcripts, sequences) = parsefa::read_fa("files/Homo_sapiens.GRCh38.cdna.all.fa");
    //let sequences = 0;

    let now = Instant::now();
    
    //splitter::contig_index_vec(sequences, contig_length);

    println!("Transcripts: {}", transcripts.len());

    let mut encoded = std::fs::read("binary.idx").unwrap();
    let mut equivalence_classes: HashMap<u32, Vec<i32>> = bincode::deserialize(&encoded[..]).unwrap();
    std::thread::spawn(move || drop(encoded));

    println!("Index size {}", equivalence_classes.len());
    println!("Load index: {}m {}s", now.elapsed().as_millis()/60000, (now.elapsed().as_millis()%60000)/1000);

    let now = Instant::now();
    let (line_count, transcript_counts) = align::read_fastq("files/SRR3534129.fastq", contig_length, &equivalence_classes, &(contig_offset as u32));
    println!("Alignment: {}m {}s", now.elapsed().as_millis()/60000, (now.elapsed().as_millis()%60000)/1000);

    std::thread::spawn(move || drop(equivalence_classes));

    let mut total_count = 0;
    let mut output = File::create("output/counts.tsv")?;
    for i in 0..transcripts.len() {
        if transcript_counts.contains_key(&(i as u32)) {
            total_count = total_count + transcript_counts.get(&(i as u32)).unwrap();
            write!(output, "{}\t{}\n", transcripts[i], transcript_counts.get(&(i as u32)).unwrap());
        }
        else{
            write!(output, "{}\t{}\n", transcripts[i], 0);
        }
    }

    println!("Total reads: {}", total_count);

    Ok(())
}

