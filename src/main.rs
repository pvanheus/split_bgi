extern crate flate2;
extern crate bio;
extern crate clap;

use bio::io::fastq::{Reader, Record, Writer};
use std::fs::File;
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use std::io;
use std::io::{BufReader, BufWriter, Error};
use clap::App;

fn split_bgi(fastq_input_path: &str) -> Result<(), Error> {
    let fastq_file = BufReader::new(File::open(fastq_input_path)?);
    let gzipped_reader = GzDecoder::new(fastq_file)?;
    let reader = Reader::new(gzipped_reader);

    let path_parts: Vec<&str> = fastq_input_path.split(".fastq.gz").collect();
    let read_1_output_path = path_parts[0].to_owned() + "_1.fastq.gz";
    let read_2_output_path = path_parts[0].to_owned() + "_2.fastq.gz";
    let mut fastq_output_1 = Writer::new(GzEncoder::new(BufWriter::new(File::create(read_1_output_path)?), Compression::default()));
    let mut fastq_output_2 = Writer::new(GzEncoder::new(BufWriter::new(File::create(read_2_output_path)?), Compression::default()));

    let iter = reader.records();
//    let mut read_2_count = 0;
    let mut read_1_count = 0;
    for result in iter {
        let record = result.unwrap();
        let description = record.desc().unwrap();
        if description.ends_with('2') {
//            read_2_count += 1;
            if let Err(e) = fastq_output_2.write_record(&record) { return Err(e) }
        } else {
            read_1_count += 1;
            let id = record.id();
            let id_parts: Vec<&str> = id.split('.').collect();
            let new_id = id_parts[0].to_owned() + "." + &read_1_count.to_string();
            let new_description = read_1_count.to_string() + "/1";
            let new_record = Record::with_attrs(
                new_id.as_str(),
                Some(new_description.as_str()),
                record.seq(),
                record.qual()
            );
            if let Err(e) = fastq_output_1.write_record(&new_record) { return Err(e) }
        }

    }
    Ok(())
}

fn main()  -> io::Result<()> {
    let matches = App::new("split_bgi")
        .version("0.1")
        .about("Splits BGISEQ-500 files from Chen et al paper")
        .author("Peter van Heusden <pvh@sanbi.ac.za>")
        .args_from_usage(
            "<INPUT>            '.fastq.gz input file with read 2 followed by read 1"
        )
        .get_matches();

    let fastq_input_path = matches.value_of("INPUT").unwrap();
    split_bgi(fastq_input_path)
}
