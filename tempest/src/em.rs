use std::collections::HashMap;
use std::time::{Instant};
use std::collections::BinaryHeap;
use std::hash::{Hash, Hasher};
use std::collections::hash_map::DefaultHasher;
use std::cmp::max;

extern crate libc;
use std::mem;

use fxhash::FxHasher;
use nohash_hasher::IntMap;
use cityhash::cityhash_1_1_1::city_hash_64;



pub fn expection_maximization2(input: Vec<Vec<u32>>, transcript_length: usize, transcript_lengths: Vec<usize>) -> Vec<f32> {
    let now = Instant::now();
    let mut ec_counts: IntMap<u32, f32> = IntMap::default();
    let mut ec_map: IntMap<u32, Vec<u32>> = IntMap::default();

    for v in input.clone() {
        let temp = listhash(v.clone());
        *ec_counts.entry(temp).or_insert(0.0) += 1.0;
        ec_map.insert(temp, v);
    }

    let mut final_counts = vec![1.0; transcript_length];

    let ec_keys: Vec<u32> = ec_counts.keys().cloned().collect::<Vec<u32>>();
    //println!("Number keys: {}", ec_keys.len());

    let mut ec_count_vec = vec![];
    let mut ec_map_vec = vec![];
    for k in 0..ec_keys.len() {
        let key = &(ec_keys[k] as u32);
        let v = ec_map.get(key).unwrap();
        
        if v.len() < 2000 {
            let mut vu = vec![];
            for i in v {
                vu.push(*i as usize);
            }
            let count = ec_counts.get(key).unwrap();
            ec_count_vec.push(count);
            ec_map_vec.push(vu);
        }
    }

    let mut new_final_counts: Vec<f32> = vec![0.0; transcript_length];
        
    for k in 0..ec_map_vec.len() {
        
        let v = &ec_map_vec[k];
        let ec_count = ec_count_vec[k];

        if ec_count > &0.0 {
            if v.len() == 1 {
                new_final_counts[v[0]] += ec_count;
            }
            else{
                for j in 0..v.len() {
                    new_final_counts[v[j]] += ec_count/v.len() as f32;
                }
            }
        }
    }
    final_counts = new_final_counts;
    
    println!("EM time: {}", now.elapsed().as_millis());

    let mut final_count_total = 0.0;
    for j in 0..final_counts.len() {
        final_count_total += final_counts[j];
    }
    //println!("total counts: {}", final_count_total);

    return final_counts;
}


pub fn expection_maximization(input: Vec<Vec<u32>>, transcript_length: usize, transcript_lengths: Vec<usize>) -> Vec<f32> {

    let mut ec_counts: IntMap<u32, f32> = IntMap::default();
    let mut ec_map: IntMap<u32, Vec<u32>> = IntMap::default();

    for v in input.clone() {
        let temp = listhash(v.clone());
        *ec_counts.entry(temp).or_insert(0.0) += 1.0;
        ec_map.insert(temp, v);
    }

    let mut final_counts = vec![1.0; transcript_length];

    let ec_keys: Vec<u32> = ec_counts.keys().cloned().collect::<Vec<u32>>();
    //println!("Number keys: {}", ec_keys.len());

    let mut ec_count_vec = vec![];
    let mut ec_map_vec = vec![];
    for k in 0..ec_keys.len() {
        let key = &(ec_keys[k] as u32);
        let v = ec_map.get(key).unwrap();
        
        if v.len() < 5000 {
            let mut vu = vec![];
            for i in v {
                vu.push(*i as usize);
            }
            let count = ec_counts.get(key).unwrap();
            ec_count_vec.push(count);
            ec_map_vec.push(vu);
        }
    }

    //println!("em classes: {}",ec_map_vec.len());

    let now = Instant::now();
    
    for i in 0..100 {
        
        let mut new_final_counts: Vec<f32> = vec![0.0; transcript_length];
        
        for k in 0..ec_map_vec.len() {
            
            let v = &ec_map_vec[k];
            let ec_count = ec_count_vec[k];

            if ec_count > &0.0 {
                if v.len() == 1 {
                    new_final_counts[v[0]] += ec_count;
                }
                else{
                    let mut weight_sum = 0.0;
                    for j in 0..v.len() {
                        //weight_sum += final_counts[v[j]]/(transcript_lengths[v[j]] as f32);
                        weight_sum += final_counts[v[j]] / (max(3, transcript_lengths[v[j]] - 200) as f32);
                    }
                    for j in 0..v.len() {
                        //let weight = final_counts[v[j]]/(transcript_lengths[v[j]] as f32);
                        let weight = final_counts[v[j]] / (max(3, transcript_lengths[v[j]] - 200) as f32);
                        let norm_weight = weight/weight_sum;
                        if norm_weight < 0.0000000001 {
                            new_final_counts[v[j]] = 0.0;
                        }
                        else{
                            new_final_counts[v[j]] += ec_count*norm_weight;
                        }
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
    //println!("total counts: {}", final_count_total);

    return final_counts;
}

fn most_frequent(array: &Vec<u32>) -> (Vec<u32>, u32) {

    let mut map:IntMap<u32, u32> = IntMap::default();

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

fn find_keys_for_value<'a>(map: &'a IntMap<u32, u32>, value: &'a u32) -> Vec<&'a u32> {
    map.iter()
        .filter_map(|(key, &val)| if &val == value { Some(key) } else { None })
        .collect()
}


pub fn listhash_old(input: Vec<u32>) -> u32 {
    let s = format!("{:?}", input);
    let mut hasher = FxHasher::default();
    hasher.write_u32(0);
    &s.hash(&mut hasher);
    return hasher.finish() as u32;
}

pub fn listhash(input: Vec<u32>) -> u32 {
    let mut hasher = FxHasher::default();
    hasher.write_u32(0);
    for &number in &input {
        hasher.write_u32(number);
    }
    hasher.finish() as u32
}