use std::cmp;
use std::fs::File;
use std::iter::FromIterator;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::error::Error;
use rand::Rng;

pub fn read_fa(input_file: &str) -> (Vec<String>, Vec<String>) {

    let f = BufReader::new(File::open(input_file).unwrap());
    
    let mut transcripts: Vec<String> = vec![];
    let mut sequences: Vec<String> = vec![];
    let mut seq = "".to_owned();
    let mut first = 1;

    let mut counter = 0;

    for (i, line) in f.lines().enumerate() {

        let lw = line.unwrap();
        counter = counter + 1;
        if counter > 500000 {
            //break;
        }

        if &lw[..1] == ">" {
            
            if first == 0 {
                sequences.push(seq.to_string().clone());
            }

            first = 0;
            
            let v: Vec<&str> = lw.split(" ").collect();
            let ensembl_transcript_id = v[0].to_string().replace(">", "");
            transcripts.push(ensembl_transcript_id);
            seq = "".to_owned();
        }
        else {
            seq = seq.to_owned();
            seq.push_str(&lw);
        }
    }

    sequences.push(seq.to_string().clone());
    
    return (transcripts, sequences);
}
