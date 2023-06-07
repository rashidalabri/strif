use clap::{Parser, Subcommand};
use extract::extract;
use merge::merge;
use profile::profile;
use std::path::PathBuf;

use crate::utils::get_default_out_path;

pub mod extract;
pub mod merge;
pub mod profile;
pub mod utils;

#[derive(Debug, Parser)]
#[clap(author, version, about, long_about = None)]
pub struct App {
    #[clap(short, long, default_value = "2")]
    verbosity: usize,

    #[clap(subcommand)]
    command: Command,
}

#[derive(Debug, Subcommand)]
enum Command {
    /// Extracts repeat sequences from an ExpansionHunter BAMlet
    Extract {
        /// The path to the ExpansionHunter BAMlet
        bamlet: PathBuf,

        /// The path to write the repeat sequences to. Defaults to the same directory as the BAMlet.
        output: Option<PathBuf>,
    },
    /// Profiles extracted repeat sequences for interruptions
    Profile {
        /// The path to the file containing repeat sequences
        repeat_seqs: PathBuf,

        /// The path to a JSON file containing the catalog of repeat loci
        str_catalog: PathBuf,

        /// Output visual alignments. Default is false.
        #[clap(short = 'z', action)]
        visual_alignments: bool,

        /// The path to the interruption profile output file. Defaults to the same directory as the repeat sequences.
        output: Option<PathBuf>,

        /// The path to visual alignments file. Defaults to the same directory as the repeat sequences.
        output_alignments: Option<PathBuf>,

        /// Filter locus IDs using a regular expression. Defaults to None.
        /// This is useful for filtering out loci that are not of interest.
        #[clap(short = 'f', long)]
        filter: Option<String>,

        #[clap(short = 'A', default_value = "1")]
        match_score: i32,

        #[clap(short = 'B', default_value = "8")]
        mismatch_penalty: i32,

        #[clap(short = 'O', default_value = "10")]
        gap_open_penalty: i32,

        #[clap(short = 'E', default_value = "1")]
        gap_extend_penalty: i32,
    },
    /// Merges profiles from multiple BAMlets partioned by case-control status
    Merge {
        /// The path to the manifest file containing paths to BAMlets and case-control status
        manifest: PathBuf,

        /// The path to a TSV file containing the global average read depth for each sample
        /// This is needed for normalizing the interription counts.
        /// The file should have two columns: sample ID and read depth.
        /// The sample ID should match the sample ID in the manifest.
        read_depths: PathBuf,

        /// The path to the merged profile. Defaults to the same directory as the manifest.
        output: Option<PathBuf>,

        /// Filter locus IDs using a regular expression. Defaults to None.
        /// This is useful for filtering out loci that are not of interest.
        #[clap(short = 'f', long)]
        filter: Option<String>,

        /// Minimum read count to include in the merged profile. Defaults to 1.
        /// This is useful for filtering out loci with low coverage.
        #[clap(short = 'm', long, default_value = "1")]
        min_read_count: u32,

        /// The sequencing read length. Used for normalizing the interruption counts.
        #[clap(short = 'l', long, default_value = "150")]
        read_length: u32,
    },
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let app: App = App::parse();

    // Set up logging
    stderrlog::new()
        .module(module_path!())
        .verbosity(app.verbosity)
        .timestamp(stderrlog::Timestamp::Second)
        .color(stderrlog::ColorChoice::Never)
        .init()
        .unwrap();

    // Match the subcommand and call relevant function with arguments
    match app.command {
        Command::Extract { bamlet, output } => {
            let out_path: PathBuf =
                output.unwrap_or_else(|| get_default_out_path(&bamlet, "repeat_seqs", "tsv"));
            extract(bamlet, out_path)?;
        }
        Command::Profile {
            repeat_seqs,
            str_catalog,
            visual_alignments,
            output,
            output_alignments,
            match_score,
            mismatch_penalty,
            gap_open_penalty,
            gap_extend_penalty,
            filter,
        } => {
            let align_parms = utils::AlignmentScoreParams {
                match_score,
                mismatch_penalty,
                gap_open_penalty,
                gap_extend_penalty,
            };
            let out_path: PathBuf = output
                .unwrap_or_else(|| get_default_out_path(&repeat_seqs, "strif_profile", "tsv"));
            let output_alns_path: PathBuf = output_alignments
                .unwrap_or_else(|| get_default_out_path(&repeat_seqs, "viz_align", "txt"));
            profile(
                repeat_seqs,
                str_catalog,
                out_path,
                output_alns_path,
                align_parms,
                visual_alignments,
                filter,
            )?;
        }
        Command::Merge {
            manifest,
            read_depths,
            output,
            filter,
            min_read_count,
            read_length,
        } => {
            let out_path: PathBuf =
                output.unwrap_or_else(|| get_default_out_path(&manifest, "merged_profile", "tsv"));
            merge(
                manifest,
                read_depths,
                out_path,
                filter,
                min_read_count,
                read_length,
            )?;
        }
    }

    Ok(())
}
