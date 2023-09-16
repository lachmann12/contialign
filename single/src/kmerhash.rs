use std::cmp::PartialEq;
use std::ops::BitXor;
use std::ops::Shl;

pub fn hash(input: &str) -> u64 {
    let mut bits: Vec<u8> = vec![0; 64];
    let mut counter: usize = 0;
    for c in input.chars() {
        //A = 00
        //T = 01
        //C = 10
        //G = 11
        if c == 'T' {
            bits[counter+1] = 1; 
        }
        else if c == 'C' {
            bits[counter] = 1;
        }
        else if c == 'G' {
            bits[counter] = 1;
            bits[counter+1] = 1;
        }
        counter += 2;
    }

    return convert(&bits).unwrap();
}

#[derive(Debug)]
pub enum ConversionError {
    Overflow,
    NonBinaryInput,
}

pub fn convert<T: PartialEq + From<u8> + BitXor<Output=T> + Shl<Output=T> + Clone>(
    bits: &[u8],
) -> Result<T, ConversionError> {
    if bits.len() > (std::mem::size_of::<T>() * 8) {
        return Err(ConversionError::Overflow);
    }
    if bits.iter()
        .filter(|&&bit| bit != 0 && bit != 1).count() > 0 {
        return Err(ConversionError::NonBinaryInput);
    }

    Ok(bits.iter()
        .fold(T::from(0), |result, &bit| {
            (result << T::from(1)) ^ T::from(bit)
        }))
}