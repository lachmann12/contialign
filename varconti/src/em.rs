use std::collections::HashMap;
use std::time::{Instant};
use std::collections::BinaryHeap;
use std::hash::{Hash, Hasher};
use std::collections::hash_map::DefaultHasher;

extern crate libc;
use std::mem;

use fxhash::FxHasher;

pub fn expection_maximization(input: Vec<Vec<u32>>, transcript_length: usize, transcript_lengths: Vec<usize>) -> Vec<f32> {

    let mut ec_counts: HashMap<u32, f32> = HashMap::new();
    let mut ec_map: HashMap<u32, Vec<u32>> = HashMap::new();

    for v in input.clone() {
        let temp = listhash(&v);
        *ec_counts.entry(temp).or_insert(0.0) += 1.0;
        ec_map.insert(temp, v);
    }

    let mut final_counts = vec![1.0; transcript_length];

    let ec_keys: Vec<u32> = ec_counts.keys().cloned().collect::<Vec<u32>>();
    println!("Number keys: {}", ec_keys.len());

    let mut ec_count_vec = vec![];
    let mut ec_map_vec = vec![];
    for k in 0..ec_keys.len() {
        let key = &(ec_keys[k] as u32);
        let v = ec_map.get(key).unwrap();
        
        if v.len() < 2 {
            let mut vu = vec![];
            for i in v {
                vu.push(*i as usize);
            }
            let count = ec_counts.get(key).unwrap();
            ec_count_vec.push(count);
            ec_map_vec.push(vu);
        }
    }

    println!("em classes: {}",ec_map_vec.len());

    let now = Instant::now();

    for i in 0..100 {
        
        let mut new_final_counts: Vec<f32> = vec![0.0; transcript_length];
        
        for k in 0..ec_map_vec.len() {
            
            let v = &ec_map_vec[k];
            let temp_count = ec_count_vec[k];

            if v.len() == 1 {
                new_final_counts[v[0]] += temp_count;
            }
            else {
                let mut gene_sum = 0.0;
                for j in 0..v.len() {
                    gene_sum += final_counts[v[j]]/(transcript_lengths[v[j]] as f32);
                }
                
                for j in 0..v.len() {
                    let temp = (final_counts[v[j]]/gene_sum)/(transcript_lengths[v[j]] as f32);
                    if temp < 0.00000000001 {
                        new_final_counts[v[j]] = 0.0;
                    }
                    else{
                        new_final_counts[v[j]] += temp_count*temp;
                    }
                }
            }
        }

        final_counts = new_final_counts;
    }

    println!("EM time: {}", now.elapsed().as_millis());

    let mut final_count_total = 0.0;
    for j in 0..final_counts.len() {
        final_count_total += final_counts[j];
    }
    println!("total counts: {}", final_count_total);

    return final_counts;
}

fn most_frequent(array: &Vec<u32>) -> (Vec<u32>, u32) {

    let mut map:HashMap<u32, u32> = HashMap::new();

    for value in array {
        *map.entry(*value).or_insert(0) += 1;
    }

    let heap: BinaryHeap<_> = map.values().collect();
    let max = heap.peek().unwrap();

    let top_elements = find_keys_for_value(&map, max);

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

pub fn listhash(input: &Vec<u32>) -> u32 {
    let mut h = FxHasher::default();
    h.write_u32(0);
    Hash::hash_slice(input, &mut h);
    return h.finish() as u32;
}