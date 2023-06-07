use lazy_static::lazy_static;
use std::{
    fs::File,
    path::{Path, PathBuf},
};

use anyhow::{Ok, Result};
use log::{info, warn};
use regex::Regex;
use rust_htslib::{
    bam,
    bam::{record::Aux, Read},
};
use std::io::prelude::*;

pub fn extract(bamlet: PathBuf, out_path: PathBuf) -> Result<()> {
    info!("Extracting repeat sequences from BAMlet...");
    let mut out_file: File = File::create(out_path)?;
    extract_repeat_seqs(&bamlet, &mut out_file)?;
    Ok(())
}

pub fn extract_repeat_seqs(bamlet: &Path, out_file: &mut File) -> Result<()> {
    // the node id of the right flank of the repeat locus (simple repeats are 2)
    let right_flank_node_id = 2;

    // captures the auxiliary tag for the repeat locus id and the cigar strings
    // for the left flank and repeat
    let formatted_regex: String = format!(
        r"^(?P<locus_id>\w+),\d+,0\[(?P<flank>(?:\d+[MIDNSHPX=])+)\](?P<repeat>(?:\d+\[(?:\d+[MIDNSHPX=])+\])+){}\[(?:\d+[MIDNSHPX=])+\]$",
        right_flank_node_id
    );
    let re_parse_tag: Regex = Regex::new(&formatted_regex).unwrap();

    let mut bam = bam::Reader::from_path(bamlet).unwrap();

    for (i, record) in bam.records().enumerate() {
        let record = record.unwrap();
        let tag: Aux = record.aux(b"XG")?;

        let tag_str = if let Aux::String(tag_str) = tag {
            tag_str
        } else {
            warn!("Auxiliary tag for read {} is not a string, skipping...", i);
            continue;
        };

        let parsed_tag = {
            if let Some(parsed_tag) = re_parse_tag.captures(tag_str) {
                parsed_tag
            } else {
                continue;
            }
        };

        let locus_id: &str = parsed_tag.name("locus_id").unwrap().as_str();
        let left_flank_cigar = parsed_tag.name("flank").unwrap().as_str();
        let repeat_cigar = parsed_tag.name("repeat").unwrap().as_str();

        // start is equal to the sum of the operation counts in the left flank
        let repeat_start = sum_operation_counts(left_flank_cigar) as usize;
        // stop is equal to the start plus the sum of operation counts in the repeat
        let repeat_stop = repeat_start + sum_operation_counts(repeat_cigar) as usize;

        let seq_raw = record.seq().as_bytes();
        let repeat_seq = std::str::from_utf8(&seq_raw[repeat_start..repeat_stop]).unwrap();

        writeln!(out_file, "{}\t{}", locus_id, repeat_seq)?;
    }

    Ok(())
}

fn sum_operation_counts(cigar: &str) -> u32 {
    // captures the numbers associated with operations that consume the read sequence
    lazy_static! {
        static ref RE_PARSE_OPS: Regex = Regex::new(r"(\d+)(?:[MISX=])").unwrap();
    }

    RE_PARSE_OPS
        .captures_iter(cigar)
        .map(|n| n[1].to_string().parse::<u32>().unwrap())
        .sum()
}
