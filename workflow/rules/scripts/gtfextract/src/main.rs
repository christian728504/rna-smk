use clap::{value_parser, Arg, Command};
use std::path::PathBuf;
use std::io::{BufRead, BufReader, Write};
use std::fs::File;
use std::collections::{HashMap, BTreeSet};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let matches = Command::new("gtfextract")
        .about("Complementy utility for hisat2-build")
        .about("Takes a gencode GTF file as input and outputs an tab-delimited exons or splice sites file acceptable to hisat2-build")
        .subcommand_required(true)
        .subcommand(
            Command::new("splice-sites")
            .about("Extract splice sites")
            .arg(
                Arg::new("gtf")
                    .long("gtf")
                    .help("Input gencode GTF file")
                    .value_parser(value_parser!(PathBuf))
                    .required(true)
            )
            .arg(
                Arg::new("out")
                    .long("out")
                    .help("Output path for splice sites")
                    .value_parser(value_parser!(PathBuf))
            )
        )
        .subcommand(
            Command::new("exons")
            .about("Extract exons")
            .arg(
                Arg::new("gtf")
                    .long("gtf")
                    .help("Input gencode GTF file")
                    .value_parser(value_parser!(PathBuf))
                    .required(true)
            )
            .arg(
                Arg::new("out")
                    .long("out")
                    .help("Output path for exons file")
                    .value_parser(value_parser!(PathBuf))
            )
        )
        .get_matches();
    
    let mut extract_ss: bool = false;
    let mut extract_exons: bool = false;
    let gtf;
    let output;

    match matches.subcommand() {
        Some(("splice-sites", args)) => {
            extract_ss = true;
            gtf = args.get_one::<PathBuf>("gtf").unwrap();
            output = args.get_one::<PathBuf>("out");
        },
        Some(("exons", args)) => {
            extract_exons = true;
            gtf = args.get_one::<PathBuf>("gtf").unwrap();
            output = args.get_one::<PathBuf>("out");
        },
        _ => unreachable!("Subcommand not found")
    }
    
    let gtf_file = File::open(gtf).unwrap();
    let buf_reader = BufReader::new(gtf_file);
    let gtf_lines = buf_reader.lines();

    struct GenomicRegion(String, char, Vec<(u32, u32)>);

    #[derive(Eq, PartialEq, Debug, Ord, PartialOrd, Clone)]
    struct Junction(String, u32, u32, char);

    let mut trans: HashMap<String, GenomicRegion> = HashMap::new();
    let mut genes: HashMap<String, Vec<String>> = HashMap::new();

    for line in gtf_lines {
        let line = line.unwrap();
        let mut line = line.trim();

        if line.is_empty() || line.starts_with('#') {
            continue
        }

        if line.contains('#') {
            line = line.split('#').next().unwrap().trim();
        }

        let parts: Vec<&str> = line.split("\t").collect();
        let [chrom, _, feature, left, right, _, strand, _, values] = &parts[..] else { continue };
        let left: u32 = left.parse().unwrap();
        let right: u32 = right.parse().unwrap();

        if feature != &"exon" || left >= right {
            continue
        };

        let mut values_map: HashMap<String, String> = HashMap::new();
        
        for attribute in values.split(";") {
            if !attribute.is_empty() {
                let trimmed = attribute.trim();
                let (key, value) = trimmed.split_once(' ').unwrap();
                values_map.insert(key.to_string(), value.trim_matches('"').to_string());
            }
        }

        if !values_map.contains_key("gene_id") || !values_map.contains_key("transcript_id") {
            continue
        }

        let transcript_id = values_map.get("transcript_id").unwrap().to_owned();
        trans.entry(transcript_id.clone())
            .and_modify(|region| region.2.push((left, right)))
            .or_insert_with(|| {
                let interval = vec![(left, right)];
                GenomicRegion(chrom.to_string(), strand.parse().unwrap(), interval)
            });
        genes.entry(values_map.get("gene_id").unwrap().clone())
            .or_default()
            .push(transcript_id.clone())
    }

    let mut sorted_trans: HashMap<String, GenomicRegion>  = HashMap::new();
    
    // If the previous range and current range overlap by <= 5 bp, then we continue, else append
    for (transcript_id, GenomicRegion(chrom, strand, mut exons)) in trans.into_iter() {
        exons.sort();
        let mut temp_exons = vec![exons[0]];
        for i in 1..exons.len() {
            if exons.get(i).unwrap().0 as i32 - temp_exons.last().unwrap().1 as i32 <= 5 {
                temp_exons.last_mut().unwrap().1 = exons.get(i).unwrap().1
            } else {
                temp_exons.push(*exons.get(i).unwrap())
            }
        }
        sorted_trans.insert(transcript_id, GenomicRegion(chrom, strand, temp_exons));
    }

    let mut junctions: BTreeSet<Junction> = BTreeSet::new();
    if extract_ss {
        for (_, GenomicRegion(chrom, strand, exons)) in sorted_trans.into_iter() {
            for i in 1..exons.len() {
                junctions.insert(Junction(chrom.clone(), exons.get(i-1).unwrap().1, exons.get(i).unwrap().0, strand));
            }
        }
    } else if extract_exons {
        let mut temp_junctions: BTreeSet<Junction> = BTreeSet::new();
        for (_, GenomicRegion(chrom, strand, exons)) in sorted_trans.into_iter() {
            for (start, end) in exons {
                temp_junctions.insert(Junction(chrom.clone(), start, end, strand));
            }
        }
        junctions.insert(temp_junctions.pop_first().unwrap());
        for junction in temp_junctions.into_iter().collect::<Vec<Junction>>() {
            let prev_junction = junctions.last().unwrap().clone();
            if junction.0 != prev_junction.0 {
                junctions.insert(junction);
                continue
            }
            assert!(prev_junction.1 <= junction.1);
            if prev_junction.2 < junction.1 {
                junctions.insert(junction);
                continue;
            }
            if prev_junction.2 < junction.2 {
                let mut strand = prev_junction.3;
                if strand != '+' && strand != '-' {
                    strand = junction.3;
                } 
                junctions.pop_last();
                junctions.insert(Junction(prev_junction.0, prev_junction.1, junction.2, strand));
            }
        }
    }

    if let Some(output) = output {
        let mut output_file = File::create(output).unwrap();
        for Junction(chrom, mut start, mut end, strand) in junctions {
            start -= 1;
            end -= 1;
            let line = format!("{chrom}\t{start}\t{end}\t{strand}\n"); 
            output_file.write_all(line.as_bytes()).unwrap();
        } 
    } else {
        for Junction(chrom, mut start, mut end, strand) in junctions {
            start -= 1;
            end -= 1;
            println!("{chrom}\t{start}\t{end}\t{strand}"); 
        }
    }

    Ok(())
}
