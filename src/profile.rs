use std::io::prelude::*;
use std::{collections::HashMap, fs::File, path::PathBuf};

use anyhow::{Ok, Result};
use bio::alignment::pairwise::Aligner;
use bio::alignment::{Alignment, AlignmentOperation};
use log::{debug, info};

use crate::utils::AlignmentScoreParams;

struct Profile {
    interruption_counts: HashMap<String, HashMap<(String, u32), u32>>,
    read_counts: HashMap<String, u32>,
}

impl Profile {
    pub fn new() -> Self {
        Self {
            interruption_counts: HashMap::new(),
            read_counts: HashMap::new(),
        }
    }

    pub fn increment_interruption(&mut self, locus_id: &str, interruption: &str, repeat_len: u32) {
        self.interruption_counts
            .entry(locus_id.to_string())
            .or_insert_with(HashMap::new)
            .entry((interruption.to_string(), repeat_len))
            .and_modify(|count| *count += 1)
            .or_insert(1);
    }

    pub fn increment_read_count(&mut self, locus_id: &str) {
        self.read_counts
            .entry(locus_id.to_string())
            .and_modify(|count| *count += 1)
            .or_insert(1);
    }

    pub fn write_to(
        &self,
        out: PathBuf,
        motifs: &HashMap<String, String>,
        reference_regions: &HashMap<String, String>,
    ) -> Result<()> {
        let mut out_file: File = File::create(out)?;
        writeln!(
            out_file,
            "locus_id\treference_region\tmotif\tread_count\tinterruption_counts"
        )?;

        let default_read_count: u32 = 0;
        let default_interruptions: HashMap<(String, u32), u32> = HashMap::new();

        for (locus_id, motif) in motifs {
            let reference_region = reference_regions.get(locus_id).unwrap();
            let read_count = self
                .read_counts
                .get(locus_id)
                .unwrap_or(&default_read_count);
            let interruptions = self
                .interruption_counts
                .get(locus_id)
                .unwrap_or(&default_interruptions);
            let interruptions_str: String = interruptions
                .iter()
                .map(|((interruption, repeat_len), count)| {
                    format!("{}:{}:{}", interruption, repeat_len, count)
                })
                .collect::<Vec<String>>()
                .join(",");
            writeln!(
                out_file,
                "{}\t{}\t{}\t{}\t{}",
                locus_id, reference_region, motif, read_count, interruptions_str
            )?;
        }
        Ok(())
    }
}

pub fn profile(
    repeat_seqs: PathBuf,
    str_catalog: PathBuf,
    out: PathBuf,
    out_alignments: PathBuf,
    align_params: AlignmentScoreParams,
    write_alignments: bool,
    filter: Option<String>,
) -> Result<()> {
    let repeat_seqs_file: File = File::open(repeat_seqs)?;
    let mut alignments_file: Option<File> = if write_alignments {
        Some(File::create(out_alignments)?)
    } else {
        None
    };

    info!("Loading STR catalog...");
    let (motifs, reference_regions) = load_str_catalog(str_catalog, filter)?;

    let mut profile: Profile = Profile::new();

    let repeat_seqs = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_reader(repeat_seqs_file)
        .into_records();

    info!("Profiling interruptions...");

    let match_fn = |a: u8, b: u8| {
        if a == b {
            align_params.match_score
        } else {
            -align_params.mismatch_penalty
        }
    };

    let mut aligner = Aligner::new(
        -align_params.gap_open_penalty,
        -align_params.gap_extend_penalty,
        &match_fn,
    );

    for record in repeat_seqs {
        let record: csv::StringRecord = record?;
        let locus_id: String = record.get(0).unwrap().to_string();

        // skip if locus is not in STR catalog
        if !motifs.contains_key(locus_id.as_str()) {
            debug!("Skipping locus {}...", locus_id);
            continue;
        }

        let repeat_seq: String = record.get(1).unwrap().to_string();
        let motif: String = motifs.get(&locus_id).unwrap().to_string();
        let motif: Vec<u8> = motif.as_bytes().to_vec();

        let observed_seq: Vec<u8> = repeat_seq.as_bytes().to_vec();
        let pure_seq = create_pure_seq(&motif, repeat_seq.len(), 4);

        let alignment: Alignment = aligner.semiglobal(&observed_seq, &pure_seq);

        // write visual alignment to file
        if alignments_file.is_some() {
            let alignments_file = alignments_file.as_mut().unwrap();
            writeln!(alignments_file, "Locus {}:", locus_id)?;
            writeln!(
                alignments_file,
                "{}",
                alignment.pretty(&observed_seq, &pure_seq, 80)
            )?;
        }

        let interruptions: Vec<String> = find_interruptions(alignment, &observed_seq);

        profile.increment_read_count(&locus_id);

        for interruption in &interruptions {
            profile.increment_interruption(&locus_id, &interruption, observed_seq.len() as u32);
        }
    }

    info!("Writing profile to output file...");
    profile.write_to(out, &motifs, &reference_regions)?;

    info!("Done!");

    Ok(())
}

fn find_interruptions(alignment: Alignment, observed: &[u8]) -> Vec<String> {
    // Given an alignment, find the interruptions in the repeat sequence
    // by looking at the path and finding consecutive insertions or substitutions
    let path = alignment.path();
    let mut interruptions: Vec<String> = Vec::new();
    let mut interruption: Vec<u8> = Vec::new();
    for step in path.iter() {
        let (observed_idx, _, op) = step;
        if *op == AlignmentOperation::Subst || *op == AlignmentOperation::Ins {
            // if *op == AlignmentOperation::Ins {
            interruption.push(observed[*observed_idx - 1]);
        } else if !interruption.is_empty() {
            interruptions.push(String::from_utf8(interruption).unwrap());
            interruption = Vec::new();
        }
    }
    interruptions
}

fn create_pure_seq(motif: &[u8], len: usize, pad: usize) -> Vec<u8> {
    // Given a motif, create a pure sequence of the motif with length
    // len and pad the end with pad copies of the motif
    let n = len / motif.len() + 1 + pad;
    motif.repeat(n)
}

fn load_str_catalog(
    str_catalog: PathBuf,
    filter: Option<String>,
) -> Result<(HashMap<String, String>, HashMap<String, String>)> {
    // create a regex filter if provided
    let filter_regex = match filter {
        Some(filter) => Some(regex::Regex::new(&filter)?),
        None => None,
    };

    let str_catalog_file: File = File::open(str_catalog)?;
    let str_catalog: Vec<HashMap<String, String>> = serde_json::from_reader(str_catalog_file)?;
    let mut motifs: HashMap<String, String> = HashMap::new();
    let mut reference_regions: HashMap<String, String> = HashMap::new();
    for mut record in str_catalog {
        let locus_id: String = record.remove("LocusId").unwrap();

        // if a filter is provided, skip if the locus id doesn't match
        if let Some(filter_regex) = &filter_regex {
            if !filter_regex.is_match(&locus_id) {
                continue;
            }
        }

        let reference_region: String = record.remove("ReferenceRegion").unwrap();
        reference_regions.insert(locus_id.clone(), reference_region);

        let structure: String = record.remove("LocusStructure").unwrap();
        let motif = structure[1..structure.len() - 2].to_string();
        motifs.insert(locus_id, motif);
    }
    Ok((motifs, reference_regions))
}
