use std::collections::HashMap;
use std::time::{Instant};
use std::collections::BinaryHeap;
use std::hash::{Hash, Hasher};
use std::collections::hash_map::DefaultHasher;

pub fn expection_maximization(input: Vec<Vec<u32>>, transcript_length: usize) -> Vec<f32> {

    
    let now = Instant::now();

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
        let count = ec_counts.get(key).unwrap();
        ec_count_vec.push(count);
        ec_map_vec.push(v);
    }

    for i in 0..100 {
        let mut new_final_counts = vec![0.0; transcript_length];
        
        for k in 0..ec_map_vec.len() {
            
            let v = ec_map_vec[k];
            let temp_count = ec_count_vec[k];

            if v.len() == 1 {
                new_final_counts[v[0] as usize] += temp_count;
            }
            else {
                let mut gene_sum = 0.0;
                for j in 0..v.len() {
                    gene_sum += final_counts[v[j] as usize];
                }
                for j in 0..v.len() {
                    let temp = final_counts[v[j] as usize]/gene_sum;
                    if temp < 0.0000000001 {
                        new_final_counts[v[j] as usize] = 0.0;
                    }
                    else{
                        new_final_counts[v[j] as usize] += temp_count*final_counts[v[j] as usize]/gene_sum;
                    }
                }
            }
        }
        final_counts = new_final_counts;
    }

    println!("EM time: {}", now.elapsed().as_millis());

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
    let mut h = DefaultHasher::new();
    Hash::hash_slice(input, &mut h);
    return h.finish() as u32;
}