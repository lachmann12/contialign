use std::collections::HashMap;
use std::time::{Instant};
use std::collections::BinaryHeap;


pub fn expection_maximization(input: Vec<Vec<u32>>) -> Vec<f32> {

    let mut final_counts = vec![1.0; input.len()];
    
    let now = Instant::now();

    for i in 0..100 {
        let mut new_final_counts = vec![0.0; input.len()];
        for v in input.clone() {
            let mut gene_sum = 0.0;
            for j in v.clone() {
                gene_sum += final_counts[j as usize];
            }
            for j in v {
                new_final_counts[j as usize] += final_counts[j as usize]/gene_sum;
            }
        }
        final_counts = new_final_counts.clone();
    }

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