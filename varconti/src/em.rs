use std::collections::HashMap;
use std::time::{Instant};

pub fn expection_maximization(input: &HashMap<u32, Vec<u32>>) -> Vec<f32> {

    let mut test: HashMap<Vec<u32>> = HashMap::new();


    let mut input: Vec<Vec<u32>> = vec![];
    let mut temp1 = vec![];
    temp1.push(0);
    temp1.push(1);
    temp1.push(2);

    let mut temp2 = vec![];
    temp2.push(0);
    temp2.push(1);

    let mut temp3 = vec![];
    temp3.push(2);
    temp3.push(3);
    temp3.push(4);

    let mut temp4 = vec![];
    temp4.push(0);
    temp4.push(3);
    temp4.push(4);

    input.push(temp1);
    input.push(temp2);
    input.push(temp3);
    input.push(temp4);

    let mut final_counts = vec![1.0; 5];
    
    let now = Instant::now();

    for i in 0..100 {
        let mut new_final_counts = vec![0.0; 5];
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

    println!("EM time: {}", now.elapsed().as_millis());

    return final_counts;
}