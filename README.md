# strif

[![Crates.io](https://img.shields.io/crates/v/strif.svg)](https://crates.io/crates/strif)
[![CI](https://github.com/rashidalabri/strif/workflows/CI/badge.svg)](https://github.com/rashidalabri/strif/actions)

## Installation

### Download binaries

Binaries for the tool can be found under the "[Releases](https://github.com/rashidalabri/strif/releases)" tab.

### Cargo

* Install the rust toolchain in order to have cargo installed by following
  [this](https://www.rust-lang.org/tools/install) guide.
* run `cargo install strif`



## Usage

### Sequence-graph alignment

To generate a sequence-graph alignment of your sample to STR loci, use [ExpansionHunter](https://github.com/Illumina/ExpansionHunter). The tool will produce a `.realigned.bam` file for each sample. Instructions for running ExpansionHunter can be found [here](https://github.com/Illumina/ExpansionHunter/blob/master/docs/03_Usage.md). 

### Extracting repeat sequences

To extract repeat sequences from an [ExpansionHunter](https://github.com/Illumina/ExpansionHunter) BAMlet (`.realigned.bam`), run the following command. If the output is not specified, the output will be saved in the same directory as the BAMlet with a `.repeat_seqs.tsv` suffix.

```
strif extract <BAMLET> [OUTPUT] 
```

### Profiling STR interruptions

To profile STR interruptions from extracted repeat sequences (the output of `strif extract`), run the following command. The STR catalog needs to be in the same format as [these catalogs](https://github.com/Illumina/RepeatCatalogs). If the output path is not specified, the output will be saved in the same directory as the repeat sequences file with a `.strif_profile.tsv` suffix.
```
strif profile [OPTIONS] <REPEAT_SEQS> <STR_CATALOG> [OUTPUT] [OUTPUT_ALIGNMENTS]
```
#### Options
```
  -z                           Output visual alignments. Default is false
  -f, --filter <FILTER>        Filter locus IDs using a regular expression. Defaults to None. This is useful for filtering out loci that are not of interest
  -A <MATCH_SCORE>             [default: 1]
  -B <MISMATCH_PENALTY>        [default: 8]
  -O <GAP_OPEN_PENALTY>        [default: 10]
  -E <GAP_EXTEND_PENALTY>      [default: 1]
```

### Merging STR interruption profiles

To merge STR interruption profiles from multiple samples, run the following command. If the output path is not specified, the output will be saved in the same directory as the manifest file with a `.merged_profiles.tsv` suffix.
```
strif merge [OPTIONS] <MANIFEST> <READ_DEPTHS> [OUTPUT]
```

- Manifest
  - Tab-separated file with the following columns:
    - Sample ID, sample status (case or control), path to STRIF profile
  - Do not include a header
  - Example
    - ```
      DO45195_case	case	output/DO45195_case.strif_profile.tsv
      DO45195_control	control	output/DO45195_control.strif_profile.tsv
      DO45231_case	case	output/DO45231_case.strif_profile.tsv
      DO45231_control	control	output/DO45231_control.strif_profile.tsv
      ```
- Read depths
  - Tab-separated file with the following columns:
    - Sample ID, read depth
  - Do not include a header
  - Example
    - ```
      DO219580_case	73.15
      DO219580_control	34.47
      DO22836_case	69.76
      DO22836_control	35.62
      ```
#### Options
```
  -f, --filter <FILTER>
          Filter locus IDs using a regular expression. Defaults to None. This is useful for filtering out loci that are not of interest
  -m, --min-read-count <MIN_READ_COUNT>
          Minimum read count to include in the merged profile. Defaults to 1. This is useful for filtering out loci with low coverage [default: 1]
  -l, --read-length <READ_LENGTH>
          The sequencing read length. Used for normalizing the interruption counts [default: 150]
  -h, --help
```

### Prioritizing interruptions

To find interruptions that display a significant difference between case and control samples, you can use `prioritize.py` in the `scripts` directory.

The prioritization script expects Sample IDs to be formatted as follows: `<INDIVIDUAL>_<case/control>`. If a paired test is run using the `-t` option, then it is expected that each individual has exactly one case and one control file.

```
python prioritize.py <merged_profile> <output_file> <sig_output_file>
```

- Output file
  - File containing information about all tested interruption, including p-values and effect sizes
  - Does not include interruption counts
- Sig(nificant) output file
  - File containining information about all interruptions with a p-value below the cut-off
  - Includes interruption counts (helpful for plotting data)

> Note: Currently, the script does not perform multiple hypothesis test correction. It is strongly recommended to independently perform this step.

#### Options
```
  -n MIN_SAMPLES, --min-samples MIN_SAMPLES
                        Minimum number of samples per group (case or control)
  -p P_VALUE_CUTOFF, --p-value-cutoff P_VALUE_CUTOFF
                        P-value cutoff
  -t, --paired-test     Enable paired test
  -c CHUNK_SIZE, --chunk-size CHUNK_SIZE
                        Chunk size for reading merged profile
  --no-progress         Disable progress bars
```

### Generating validation datasets
You can generate simulate repeat sequences to validate and test STRIF using `generate_validation_sets.py` in the `scripts` directory. The only argument is a path to a directory, such as `datasets/` where the generated datasets will be created.

```
python generate_validation_sets.py <DATASET_DIR>
```

- Generated datasets
  - `simple`
    - Small dataset helpful for debugging
  - `no_interruption`
    - Repeat sequences containing no interruptions
  - `basic_<1-6>`
    - Small dataset useful for development
  - `comprehensive_<test, train, valid>`
    - Comprehensive dataset useful for optimizing parameters, validating and testing
  - `disjoint_<1-6>`
    - Dataset of disjoint interruptions where the interruption sequence does not include any bases from the repeat sequence
  - `intersect_<1-6>`
    - Dataset of intersecting interruptions where the interruption sequence includes at least one base from the repeat sequence
  - `insert_<1-6>`
    - Dataset of interruptions that have been inserted into the repeat sequence
  - `substitute_<1-6>`
    - Dataset of interruptions that have substituted one or more repeat sequence bases

### Calculating performance metrics
You can calculate metrics on the generated datasets using `metrics.py` in the `scripts` directory. The only argument is a path to a directory, such as `datasets/` where the generated datasets was created.

```
python metrics.py <DATASET_DIR>
```

The script will output a file `overall_stats.tsv` in the dataset directory containing a summary of metrics on each dataset.

### Optimizing alignment parameters
You can find optimal aligning parameters for `strif profile` by running `optimize.py` in the `scripts` directory. The only argument is a path to a dataset. This will be any directory within the datasets directory. It is recommended to run this on `datasets/comprehensive_train`.

```
python optimize.py <DATASET_DIR>/<NAME_OF_DATASET>
```

## License

Licensed under either of

 * Apache License, Version 2.0
   ([LICENSE-APACHE](LICENSE-APACHE) or http://www.apache.org/licenses/LICENSE-2.0)
 * MIT license
   ([LICENSE-MIT](LICENSE-MIT) or http://opensource.org/licenses/MIT)

at your option.

## Contribution

Unless you explicitly state otherwise, any contribution intentionally submitted
for inclusion in the work by you, as defined in the Apache-2.0 license, shall be
dual licensed as above, without any additional terms or conditions.

See [CONTRIBUTING.md](CONTRIBUTING.md).
