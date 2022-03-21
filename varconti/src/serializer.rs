
use rkyv::{Archive, Deserialize, Serialize};


use rkyv::{
  ser::{serializers::{AlignedSerializer, AllocSerializer}, Serializer},
  util::AlignedVec,
};

use std::collections::HashMap;
use std::fs;
use std::io::Write;

use terminal_spinners::{SpinnerBuilder, BOUNCE};

pub fn serialize(filename: &String, transcripts: &Vec<String>, eq_classes: &HashMap<u64, u32>, eq_elements: &HashMap<u32, Vec<u32>>, transcript_kmers: &HashMap<u32, u32>, kmer_length: &u32, index_version: &u32, transcript_length: &Vec<usize>) -> std::io::Result<()> {

    let handle = SpinnerBuilder::new().spinner(&BOUNCE).text(" Serializing transcripts ... ").start();

    let mut file = fs::OpenOptions::new()
      .write(true)
      .truncate(true)
      .create(true)
      .open(filename)
      .unwrap();

    let bytes_size = rkyv::to_bytes::<_, 64>(index_version).unwrap();
    file.write_all(&bytes_size)?;

    let bytes_size = rkyv::to_bytes::<_, 64>(kmer_length).unwrap();
    file.write_all(&bytes_size)?;

    let bytes = rkyv::to_bytes::<_, 64>(transcripts).unwrap();
    let bytes_size = rkyv::to_bytes::<_, 64>(&(bytes.len() as u64)).unwrap();
    file.write_all(&bytes_size)?;
    file.write_all(&bytes)?;

    handle.text(" Serializing EQ elements ... ");
    let bytes = rkyv::to_bytes::<_, 64>(eq_elements).unwrap();
    let bytes_size = rkyv::to_bytes::<_, 64>(&(bytes.len() as u64)).unwrap();
    file.write_all(&bytes_size)?;
    file.write_all(&bytes)?;

    handle.text(" Serializing Transcript kmer count ... ");
    let bytes = rkyv::to_bytes::<_, 64>(transcript_kmers).unwrap();
    let bytes_size = rkyv::to_bytes::<_, 64>(&(bytes.len() as u64)).unwrap();
    file.write_all(&bytes_size)?;
    file.write_all(&bytes)?;

    // test serialize
    handle.text(" testing Serializing EQ classes ... ");
    let mut serializer = AllocSerializer::<256>::default();
    serializer.serialize_value(eq_classes).unwrap();
    let bytes = serializer.into_serializer().into_inner();

    handle.text(" Serializing EQ classes ... ");
    let bytes = rkyv::to_bytes::<_, 64>(eq_classes).unwrap();
    let bytes_size = rkyv::to_bytes::<_, 64>(&(bytes.len() as u64)).unwrap();
    file.write_all(&bytes_size)?;
    file.write_all(&bytes)?;

    handle.text(" Serializing transcript length ... ");
    let bytes = rkyv::to_bytes::<_, 64>(transcript_length).unwrap();
    let bytes_size = rkyv::to_bytes::<_, 64>(&(bytes.len() as u64)).unwrap();
    file.write_all(&bytes_size)?;
    file.write_all(&bytes)?;

    handle.text(" Index serialized ");
    handle.stop_and_clear();

    Ok(())
}

pub fn deserialize(filename: &String) -> (u32, Vec<String>,  HashMap<u32, Vec<u32>>, HashMap<u64, u32>, HashMap<u32, u32>, Vec<usize>) {
    let handle = SpinnerBuilder::new().spinner(&BOUNCE).text(" Loading index ... ").start();

    let mbytes = std::fs::read(filename).unwrap();
    //let archived = unsafe { rkyv::archived_root::<u32>(&mbytes[0..4]) };
    //let index_version: u32 = archived.deserialize(&mut rkyv::Infallible).unwrap();

    let archived = unsafe { rkyv::archived_root::<u32>(&mbytes[4..8]) };
    let kmer_size: u32 = archived.deserialize(&mut rkyv::Infallible).unwrap();

    handle.text(" Initializing transcripts ... ");
    let archived = unsafe { rkyv::archived_root::<u64>(&mbytes[8..16]) };
    let blength: u64 = archived.deserialize(&mut rkyv::Infallible).unwrap();
    let transcript_end = (16+blength as usize) as usize;
    let archived = unsafe { rkyv::archived_root::<Vec<String>>(&mbytes[16..transcript_end]) };
    let transcripts: Vec<String> = archived.deserialize(&mut rkyv::Infallible).unwrap();
    
    handle.text(" Initializing EQ elements ... ");
    let archived = unsafe { rkyv::archived_root::<u64>(&mbytes[transcript_end..(transcript_end+8)]) };
    let blength: u64 = archived.deserialize(&mut rkyv::Infallible).unwrap();
    let elements_end = (transcript_end+blength as usize +8) as usize;
    let archived = unsafe { rkyv::archived_root::<HashMap<u32, Vec<u32>>>(&mbytes[(transcript_end+8)..elements_end]) };
    let eq_elements: HashMap<u32, Vec<u32>> = archived.deserialize(&mut rkyv::Infallible).unwrap();
    
    handle.text(" Initializing Transcript kmer counts ... ");
    let archived = unsafe { rkyv::archived_root::<u64>(&mbytes[elements_end..(elements_end+8)]) };
    let blength: u64 = archived.deserialize(&mut rkyv::Infallible).unwrap();
    let trans_kmer_end = (elements_end+blength as usize +8) as usize;
    let archived = unsafe { rkyv::archived_root::<HashMap<u32, u32>>(&mbytes[(elements_end+8)..trans_kmer_end]) };
    let transcrip_kmer_count: HashMap<u32, u32> = archived.deserialize(&mut rkyv::Infallible).unwrap();
    
    handle.text(" Initializing EQ classes ... ");
    let archived = unsafe { rkyv::archived_root::<u64>(&mbytes[trans_kmer_end..(trans_kmer_end+8)]) };
    let blength: u64 = archived.deserialize(&mut rkyv::Infallible).unwrap();
    let classes_end = (trans_kmer_end+blength as usize +8) as usize;
    let archived = unsafe { rkyv::archived_root::<HashMap<u64, u32>>(&mbytes[(trans_kmer_end+8)..classes_end]) };
    let eq_classes: HashMap<u64, u32> = archived.deserialize(&mut rkyv::Infallible).unwrap();

    handle.text(" Initializing Transcript length ... ");
    let archived = unsafe { rkyv::archived_root::<u64>(&mbytes[classes_end..(classes_end+8)]) };
    let blength: u64 = archived.deserialize(&mut rkyv::Infallible).unwrap();
    let transcript_length_end = (classes_end+blength as usize +8) as usize;
    let archived = unsafe { rkyv::archived_root::<Vec<usize>>(&mbytes[(classes_end+8)..transcript_length_end]) };
    let transcript_length: Vec<usize> = archived.deserialize(&mut rkyv::Infallible).unwrap();
    
    handle.text(" Index serialized ");
    handle.stop_and_clear();

    return (kmer_size, transcripts, eq_elements, eq_classes, transcrip_kmer_count, transcript_length);
}