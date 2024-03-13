//use rust_htslib::{bam, bam::Read};
use bam;
use bam::BamReader;
use std::collections::{HashMap};
use std::env;
use std::fs::File;
use std::io::BufWriter;
use std::io::Write;

#[derive(Clone, Copy, PartialEq, Eq, Hash)]
struct Location {
    contig_num: usize,
    position: usize,
}

fn main() {

    // Collect command line arguments
    let args: Vec<String> = env::args().collect();
    if args.len() != 3 {
        show_usage_and_exit();
    }
    let input_bam_path = &args[1];
    let output_prefix = &args[2];


    // Try to create a bam reader
    let bam_reader = bam::BamReader::from_path(input_bam_path, 1).unwrap_or_else(|error| {
        panic!("Could not read the input bam file: {:?}", error)
    });

    let contig_names = get_sequence_names_from_header(&bam_reader);

    let mapq_sums = collect_stats(bam_reader).expect("Could not collect stats from bam");

    // Calculate average mapqs for each location
    let mut average_mapqs: Vec<(Location, i32)> = Vec::new();
    for (location, (sum, count)) in &mapq_sums {
        let average_mapq = (sum / *count as f64).round() as i32;
        average_mapqs.push((*location, average_mapq));
    }

    // Sort the average_mapqs by first contig_num and then position
    average_mapqs.sort_by_key(|&(location, _)| (location.contig_num, location.position));


    // Output bedgraph files
    write_mapq_bedgraph(&average_mapqs, &contig_names, format!("{}.mapq.bedgraph", output_prefix));
    write_depth_bedgraph(&mapq_sums, &average_mapqs, &contig_names, format!("{}.dp.bedgraph", output_prefix));
}





fn show_usage_and_exit() {
    println!("USAGE: bdp BAM_FILE OUTPUT_PREFIX");
    std::process::exit(1);
}

fn collect_stats(bam_reader: BamReader<File>) -> Result<HashMap<Location, (f64, usize)>, Box<dyn std::error::Error>> {
    let mut mapq_sums: HashMap<Location, (f64, usize)> = HashMap::new();
    for r in bam_reader {
        let record = r?;

        // Skip unmapped reads
        if record.start() < 0 {
            continue
        }

        let mut curr_pos = record.start();
        for (cigar_len, cigar_op) in record.cigar().iter() {
            if cigar_op.consumes_ref() {
                for _ in 1..cigar_len + 1 {

                    // Create Location struct
                    let location = Location {
                        contig_num: record.ref_id() as usize,
                        position: curr_pos as usize
                    };

                    // Add up MAPQ scores and count reads for the location
                    let (sum, count) = mapq_sums.entry(location).or_insert((0.0, 0));
                    *sum += record.mapq() as f64;
                    *count += 1;

                    //positions_covered.insert(location);

                    curr_pos += 1;
                }
            }
        }
    }
    return Ok(mapq_sums)
}



fn get_sequence_names_from_header(bam_reader: &BamReader<File>) -> HashMap<usize, String> {
    bam_reader.header().reference_names().iter().enumerate()
        .map(|(i, name)| (i, name.to_string()))
        .collect()
}



fn write_depth_bedgraph(mapq_sums: &HashMap<Location, (f64, usize)>, average_mapqs: &Vec<(Location, i32)>, contig_names: &HashMap<usize, String>, output_path: String) {
    let dp_file = File::create(output_path).unwrap();
    let mut dp_file_writer = BufWriter::new(&dp_file);

    let mut first = true;
    let mut region_start = 0;
    let mut prev_contig = 0;
    let mut prev_position = 0;
    let mut prev_dp = 0;
    for (location, _) in average_mapqs {
        let (_, dp) = mapq_sums.get(&location).unwrap();
        if &prev_dp != dp || prev_contig != location.contig_num {
            if !first {
                writeln!(&mut dp_file_writer, "{}\t{}\t{}\t{}", contig_names.get(&prev_contig).unwrap(), region_start, prev_position+1, prev_dp).expect("Could not write to dp file");
            }
            first = false;
            region_start = location.position;
        }
        prev_dp = *dp;
        prev_contig = location.contig_num;
        prev_position = location.position
    }
    writeln!(&mut dp_file_writer, "{}\t{}\t{}\t{}", contig_names.get(&prev_contig).unwrap(), region_start, prev_position+1, prev_dp).expect("Could not write to dp file");
}



fn write_mapq_bedgraph(average_mapqs: &Vec<(Location, i32)>, contig_names: &HashMap<usize, String>, output_path: String) {
    let mapq_file = File::create(output_path).unwrap();
    let mut mapq_file_writer = BufWriter::new(&mapq_file);

    let mut region_start = 0;
    let mut prev_mapq = 0;
    let mut prev_contig = 0;
    let mut prev_position = 0;
    let mut first = true;
    for (location, mapq) in average_mapqs {
        if &prev_mapq != mapq || prev_contig != location.contig_num {
            if !first {
                writeln!(&mut mapq_file_writer, "{}\t{}\t{}\t{}", contig_names.get(&prev_contig).unwrap(), region_start, prev_position+1, prev_mapq).expect("Could not write to mapq file");
            }
            first = false;
            region_start = location.position;
        }
        prev_mapq = *mapq;
        prev_contig = location.contig_num;
        prev_position = location.position
    }
    writeln!(&mut mapq_file_writer, "{}\t{}\t{}\t{}", contig_names.get(&prev_contig).unwrap(), region_start, prev_position+1, prev_mapq).expect("Could not write to mapq file");
}
