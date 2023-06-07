import argparse
from collections import defaultdict

import pandas as pd
import scipy.stats
import numpy as np
from numpy import mean, sqrt, std
from tqdm.auto import tqdm


def cohen_d(cases, controls):
    """
    Source: https://stackoverflow.com/questions/21532471/how-to-calculate-cohens-d-in-python
    """
    nx = len(cases)
    ny = len(controls)
    dof = nx + ny - 2
    return (mean(cases) - mean(controls)) / sqrt(
        ((nx - 1) * std(cases, ddof=1) ** 2 + (ny - 1) * std(controls, ddof=1) ** 2)
        / dof
    )


def file_len(fname):
    with open(fname) as f:
        for i, _ in enumerate(f):
            pass
    return i + 1


def main():
    # Create the argument parser
    parser = argparse.ArgumentParser(description="Prioritize interruptions")
    parser.add_argument("merged_profile", type=str, help="Path to merged profile")
    parser.add_argument(
        "output_file", type=str, help="Path to output file containing all interruptions"
    )
    parser.add_argument(
        "sig_output_file",
        type=str,
        help="Path to output file containing interruptions that pass the p-value cutoff with corresponding normalized counts",
    )
    parser.add_argument(
        "-n",
        "--min-samples",
        type=int,
        default=2,
        help="Minimum number of samples per group (case or control)",
    )
    parser.add_argument(
        "-p", "--p-value-cutoff", type=float, default=0.05, help="P-value cutoff"
    )
    parser.add_argument(
        "-t", "--paired-test", action="store_true", help="Enable paired test"
    )
    parser.add_argument(
        "-c",
        "--chunk-size",
        type=int,
        default=5000,
        help="Chunk size for reading merged profile",
    )
    parser.add_argument(
        "--no-progress",
        action="store_true",
        help="Disable progress bars",
    )

    # Parse the arguments
    args = parser.parse_args()

    # Access the parsed arguments
    merged_profile_path = args.merged_profile
    output_file_path = args.output_file
    sig_output_file_path = args.sig_output_file
    min_samples = args.min_samples
    p_value_cutoff = args.p_value_cutoff
    paired_test = args.paired_test
    chunk_size = args.chunk_size
    progress_bar = args.no_progress

    # Load merged profile
    merged_profile = pd.read_csv(merged_profile_path, sep="\t", chunksize=chunk_size)
    num_chunks = int((file_len(merged_profile_path) - 1) / chunk_size + 1)

    if paired_test:
        print("Paired test enabled")
    else:
        print("Paired test disabled")

    n_skipped_interruptions = 0
    n_skipped_loci = 0

    output = []

    for chunk in tqdm(
        merged_profile,
        total=num_chunks,
        desc="Chunks",
        position=0,
        disable=progress_bar,
    ):
        # Loop over the loci in the merged profile chunk
        for _, row in tqdm(
            chunk.iterrows(),
            total=len(chunk),
            desc="Loci",
            position=1,
            leave=False,
            disable=progress_bar,
        ):
            # We exclude loci that do not have enough samples included
            # We do an initial check here to avoid unnecessary computation
            read_counts = row["read_counts"].split(",")
            if len(read_counts) < min_samples * 2:
                n_skipped_loci += 1
                continue

            # We find the donors that are included in the analysis of this locus by
            # looking at which samples have read counts
            case_donors = []
            control_donors = []

            for entry in read_counts:
                sample_id, count = entry.split(":")
                donor_id, status = sample_id.split("_")

                if status == "case":
                    if donor_id in case_donors:
                        raise Exception(f"Duplicate case sample: {sample_id}")
                    case_donors.append(donor_id)
                elif status == "control":
                    if donor_id in control_donors:
                        raise Exception(f"Duplicate control sample: {sample_id}")
                    control_donors.append(donor_id)
                else:
                    raise Exception(f"Unknown sample status: {sample_id}")

            # If we're doing a paired test, exclude donors that do not have a paired case/control sample
            if paired_test:
                paired_donors = set(case_donors).intersection(set(control_donors))
                paired_donors = list(paired_donors)
                case_donors = paired_donors
                control_donors = paired_donors

            case_donors.sort()
            control_donors.sort()

            # We skip loci that do not have enough donors included
            if len(case_donors) < min_samples or len(control_donors) < min_samples:
                n_skipped_loci += 1
                continue

            # Create a dict of interruption unit -> status -> donor_id -> count
            interruptions = defaultdict(
                lambda: {
                    "case": {d: 0 for d in case_donors},
                    "control": {d: 0 for d in control_donors},
                }
            )

            intrp_units_to_skip = set()

            # Fill the interruptions dict from the interruption_counts column
            # which contains info for all interruptions in the locus
            for entry in row["interruption_counts"].split(","):
                sample_id, intrpt_unit, count = entry.split(":")
                donor_id, status = sample_id.split("_")
                if status not in ["case", "control"]:
                    raise Exception("Unknown sample status: {}".format(sample_id))
                elif status == "case" and donor_id not in case_donors:
                    continue
                elif status == "control" and donor_id not in control_donors:
                    continue

                count = float(count)
                if (
                    pd.isna(count)
                    or count == float("inf")
                    or count == float("-inf")
                    or count < 0
                ):
                    intrp_units_to_skip.add(intrpt_unit)

                interruptions[intrpt_unit][status][donor_id] = count

            # For each interruption unit, we calculate the p-value and cohen's d
            for intrpt_unit, counts in tqdm(
                interruptions.items(),
                total=len(interruptions),
                desc="Interruptions",
                position=2,
                leave=False,
                disable=progress_bar,
            ):
                # Skip interruptions that have invalid counts
                if intrpt_unit in intrp_units_to_skip:
                    n_skipped_interruptions += 1
                    continue

                case_counts = [counts["case"][d] for d in case_donors]
                control_counts = [counts["control"][d] for d in control_donors]

                # calculate p-value
                if paired_test:
                    _, p_value = scipy.stats.wilcoxon(
                        case_counts, control_counts, alternative="two-sided"
                    )
                else:
                    _, p_value = scipy.stats.mannwhitneyu(
                        case_counts, control_counts, alternative="two-sided"
                    )

                # Calculate cohen's d
                cohen_d_value = cohen_d(case_counts, control_counts)

                if pd.isna(p_value):
                    print(
                        f"Warning: NaN p-value for {row['locus_id']} and '{intrpt_unit}' interruption. Skipping..."
                    )
                    n_skipped_interruptions += 1
                    continue
                elif pd.isna(cohen_d_value):
                    print(
                        f"Warning: NaN Cohen's d for {row['locus_id']} and '{intrpt_unit}' interruption. Skipping..."
                    )
                    n_skipped_interruptions += 1
                    continue

                # We only output the counts if the p-value is below the cutoff
                if p_value < p_value_cutoff:
                    read_counts_str = row["read_counts"]
                    interruption_counts_str = []
                    for d in case_donors:
                        interruption_counts_str.append(f"{d}_case:{counts['case'][d]}")
                    for d in control_donors:
                        interruption_counts_str.append(
                            f"{d}_control:{counts['control'][d]}"
                        )
                    interruption_counts_str = ",".join(interruption_counts_str)
                else:
                    read_counts_str = ""
                    interruption_counts_str = ""

                output.append(
                    (
                        row["locus_id"],
                        row["reference_region"],
                        row["motif"],
                        intrpt_unit,
                        len(case_donors),
                        len(control_donors),
                        p_value,
                        cohen_d_value,
                        read_counts_str,
                        interruption_counts_str,
                    )
                )

    # create a dataframe from the output
    output_cols = [
        "locus_id",
        "reference_region",
        "motif",
        "interruption",
        "n_case",
        "n_control",
        "p_value",
        "cohen_d",
        "read_counts",
        "interruption_counts",
    ]
    output_df = pd.DataFrame(output, columns=output_cols)

    # sort by p-value
    output_df.sort_values(by=["p_value"], inplace=True, ignore_index=True)

    # write to file
    output_df[output_cols[:-2]].to_csv(output_file_path, sep="\t", index=False)

    # only keep interruptions with p-value < cutoff
    output_df = output_df[output_df["p_value"] < p_value_cutoff]

    # write to file
    output_df.to_csv(sig_output_file_path, sep="\t", index=False)

    print(f"Skipped {n_skipped_loci} loci")
    print(f"Skipped {n_skipped_interruptions} interruptions")
    print("Done!")


main()
