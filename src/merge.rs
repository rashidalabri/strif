use std::io::prelude::*;
use std::{collections::HashMap, fs::File, path::PathBuf};

use anyhow::{Ok, Result};
use log::{info, warn};

type LocusId = String;
type Motif = String;
type Interruption = String;
type SampleId = String;
type ReferenceRegion = String;
type Count = u32;
type NormCount = f64;

type LocusInterruptionCounts = HashMap<(SampleId, Interruption), NormCount>;
type InterruptionCounts = HashMap<LocusId, LocusInterruptionCounts>;
type ReadCounts = HashMap<LocusId, Vec<(SampleId, Count)>>;

struct MergedProfile {
    interruption_counts: InterruptionCounts,
    read_counts: ReadCounts,
    motifs: HashMap<LocusId, Motif>,
    reference_regions: HashMap<LocusId, ReferenceRegion>,
}

impl MergedProfile {
    pub fn new() -> Self {
        Self {
            interruption_counts: HashMap::new(),
            read_counts: HashMap::new(),
            motifs: HashMap::new(),
            reference_regions: HashMap::new(),
        }
    }

    pub fn increment_interruption(
        &mut self,
        locus_id: &str,
        sample_id: &str,
        interruption: &str,
        count: f64,
    ) {
        self.interruption_counts
            .entry(locus_id.to_string())
            .or_insert_with(HashMap::new)
            .entry((sample_id.to_string(), interruption.to_string()))
            .and_modify(|c| *c += count)
            .or_insert(count);
    }

    pub fn add_read_count(&mut self, locus_id: &str, sample_id: &str, count: u32) {
        self.read_counts
            .entry(locus_id.to_string())
            .or_insert_with(Vec::new)
            .push((sample_id.to_string(), count));
    }

    pub fn add_reference_region(&mut self, locus_id: &str, reference_region: &str) {
        self.reference_regions
            .entry(locus_id.to_string())
            .or_insert(reference_region.to_string());
    }

    pub fn add_motif(&mut self, locus_id: &str, motif: &str) {
        self.motifs
            .entry(locus_id.to_string())
            .or_insert(motif.to_string());
    }

    pub fn write_to(&self, out: PathBuf) -> Result<()> {
        let mut out_file: File = File::create(out)?;
        writeln!(
            out_file,
            "locus_id\treference_region\tmotif\tread_counts\tinterruption_counts"
        )?;
        let default_interruption_counts: LocusInterruptionCounts = HashMap::new();
        for (locus_id, motif) in &self.motifs {
            let reference_region: &String = self.reference_regions.get(locus_id).unwrap();
            let read_counts: &Vec<(String, u32)> = self.read_counts.get(locus_id).unwrap();
            let interruption_counts = self
                .interruption_counts
                .get(locus_id)
                .unwrap_or(&default_interruption_counts);

            let read_counts_str = read_counts
                .iter()
                .map(|(sample_id, count)| format!("{}:{}", sample_id, count))
                .collect::<Vec<String>>()
                .join(",");
            let interruption_counts_str = interruption_counts
                .iter()
                .map(|((sample_id, interruption), count)| {
                    format!("{}:{}:{}", sample_id, interruption, count)
                })
                .collect::<Vec<String>>()
                .join(",");
            writeln!(
                out_file,
                "{}\t{}\t{}\t{}\t{}",
                locus_id, reference_region, motif, read_counts_str, interruption_counts_str
            )?;
        }
        Ok(())
    }
}

pub fn merge(
    manifest: PathBuf,
    read_depths: PathBuf,
    out_path: PathBuf,
    filter: Option<String>,
    min_read_count: u32,
    read_len: u32,
) -> Result<()> {
    info!("Merging profiles from manifest...");

    // create a regex filter if provided
    let filter_regex = match filter {
        Some(filter) => Some(regex::Regex::new(&filter)?),
        None => None,
    };

    // load manifest, which is a TSV with columns: sample, case_control, profile_path (no headers)
    let mut profiles: Vec<(SampleId, PathBuf)> = Vec::new();
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_path(manifest)?;
    for result in reader.records() {
        let record = result?;
        let sample_id = record.get(0).unwrap().to_string();
        let profile_path: PathBuf = PathBuf::from(record.get(2).unwrap());
        profiles.push((sample_id, profile_path));
    }

    // load the read depths file, which is a TSV with columns: sample, read_depth (no headers)
    let mut read_depths_map: HashMap<SampleId, f64> = HashMap::new();
    let mut read_depths_reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_path(read_depths)?;
    for result in read_depths_reader.records() {
        let record = result?;
        let sample_id = record.get(0).unwrap().to_string();
        let read_depth: f64 = record.get(1).unwrap().parse::<f64>()?;
        read_depths_map.insert(sample_id, read_depth);
    }

    // open each profile and add to merged profile
    let mut merged_profile = MergedProfile::new();
    for (sample_id, profile_path) in profiles {
        info!("Processing {} profile...", sample_id);
        let mut reader: csv::Reader<File> = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(true)
            .from_path(profile_path)?;
        for result in reader.records() {
            let record: csv::StringRecord = result?;
            let locus_id: &str = record.get(0).unwrap();

            // skip if locus_id does not match filter
            if let Some(filter_regex) = &filter_regex {
                if !filter_regex.is_match(locus_id) {
                    continue;
                }
            }

            // skip if read count is below minimum otherwise add to merged profile
            let read_count: u32 = record.get(3).unwrap().parse::<u32>()?;
            if read_count < min_read_count {
                continue;
            } else {
                merged_profile.add_read_count(locus_id, &sample_id, read_count);
            }

            // add reference region and motif to merged profile
            let reference_region: &str = record.get(1).unwrap();
            merged_profile.add_reference_region(locus_id, reference_region);
            let motif: &str = record.get(2).unwrap();
            merged_profile.add_motif(locus_id, motif);

            let interruption_counts_str: &str = record.get(4).unwrap();
            if interruption_counts_str.len() == 0 {
                continue;
            }

            let interruption_counts: Vec<&str> = interruption_counts_str.split(",").collect();
            for interruption_count in interruption_counts {
                let interruption_count: Vec<&str> = interruption_count.split(":").collect();
                let interruption: &str = interruption_count[0];
                let repeat_len: u32 = interruption_count[1].parse::<u32>()?;
                if repeat_len == 0 || repeat_len > read_len {
                    warn!("Sample {} has an invalid repeat length={} for {} with a '{}' interruption. Read length={}.", sample_id, repeat_len, locus_id, interruption, read_len);
                }
                let count: u32 = interruption_count[2].parse::<u32>()?;
                let read_depth = read_depths_map.get(&sample_id).unwrap().clone();
                let norm_count: f64 =
                    norm_interruption_count(count, read_len, repeat_len, read_depth);
                if norm_count.is_infinite() || norm_count.is_nan() || norm_count < 0.0 {
                    warn!(
                        "Sample {} has an invalid normalized count={} for {} with a '{}' interruption. Raw count={}, read length={}, repeat length={}, read depth={}.",
                        sample_id, norm_count, locus_id, interruption, count, read_len, repeat_len, read_depth
                    );
                }
                merged_profile.increment_interruption(
                    locus_id,
                    &sample_id,
                    interruption,
                    norm_count,
                );
            }
        }
    }

    merged_profile.write_to(out_path)?;

    Ok(())
}

fn norm_interruption_count(count: u32, read_len: u32, repeat_len: u32, read_depth: f64) -> f64 {
    let num_possible_start: u32 = read_len - repeat_len + 1;
    let expected_num_reads: f64 = num_possible_start as f64 * read_depth;
    (count as f64) / expected_num_reads
}
