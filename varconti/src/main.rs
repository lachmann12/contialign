use std::time::{Instant};

use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};

use std::collections::HashMap;
use std::collections::HashSet;

use std::{thread, time};

mod parsefa;
mod splitter;

use std::fs::File;
use std::io::Read;

extern crate bincode;
use std::iter::FromIterator;
use rand::Rng;

// cargo build --release
// "target/release/varconti.exe"

fn main() {

    println!("Start");

    let mut hasher = DefaultHasher::new();
    let mut contig_length = 22;

    //let (transcripts, sequences) = parsefa::read_fa("files/Homo_sapiens.GRCh38.cdna.all.fa");

    let now = Instant::now();
    
    //splitter::contig_index_vec(sequences, contig_length);
    
    //println!("Transcripts: {}", transcripts.len());

    let mut encoded = std::fs::read("binary.idx").unwrap();
    let mut equivalence_classes: HashMap<u32, Vec<i32>> = bincode::deserialize(&encoded[..]).unwrap();

    println!("Index size {}", equivalence_classes.len());
    println!("Building index: {}m {}s", now.elapsed().as_millis()/60000, (now.elapsed().as_millis()%60000)/1000);

    let testkeys = Vec::from_iter(equivalence_classes.keys());

    
    let now = Instant::now();
    let mut u = &vec![];
    
    let mut total = 0;
    for i in 1..100000000 {
        let rr = rand::thread_rng().gen_range(0..10000000);
        u = equivalence_classes.get(testkeys[rr]).unwrap();
        total = total + u.len();
    }
    println!("query time: {}ms", now.elapsed().as_millis());
    println!("result: {:?}", u);
    println!("result: {}", total);
}
