from pathlib import Path
import subprocess
import pandas as pd
import sys


def run_command(
    repeat_seqs,
    str_catalog,
    output,
    match_score=None,
    mismatch_score=None,
    gap_open_score=None,
    gap_extend_score=None,
    visualize=False,
):
    command = f"cargo run --release -- profile {repeat_seqs} {str_catalog} {output}"
    if match_score is not None:
        command += f" -A={match_score}"
    if mismatch_score is not None:
        command += f" -B={mismatch_score}"
    if gap_open_score is not None:
        command += f" -O={gap_open_score}"
    if gap_extend_score is not None:
        command += f" -E={gap_extend_score}"
    if visualize:
        command += " -z"
    subprocess.run(command, shell=True, stderr=subprocess.DEVNULL)
    return command


def create_stat_df(truth_path, profile_path):
    truth = pd.read_csv(truth_path, sep="\t")
    profile = pd.read_csv(profile_path, sep="\t")

    truth = truth.rename(columns={"interruption": "true_interruption"})
    profile = profile.rename(columns={"interruption": "pred_interruption"})

    merged = pd.merge(truth, profile, on=["locus_id", "motif"], how="outer")
    merged = merged.fillna("")
    merged = merged.drop(columns=["seq", "count", "total_count"])

    merged["exact_match"] = merged["true_interruption"] == merged["pred_interruption"]
    merged["exact_match"] = merged["exact_match"].astype(int)

    # inexact match is true if true interruption is a substring of predicted interruption
    merged["inexact_match"] = merged.apply(
        lambda row: row["true_interruption"] in row["pred_interruption"], axis=1
    )
    merged["inexact_match"] = merged["inexact_match"].astype(int)

    return merged


def compute_mean_precision(stat_df, match_col):
    mean_precision = (
        stat_df.groupby("locus_id").agg({match_col: "mean"})[match_col].mean()
    )
    return mean_precision


def compute_mean_recall(stat_df, match_col):
    mean_precision = (
        stat_df.groupby("locus_id").agg({match_col: "sum"})[match_col].mean()
    )
    return mean_precision


def compute_metrics(stats_df):
    exact_mean_precision = compute_mean_precision(stats_df, "exact_match")
    exact_mean_recall = compute_mean_recall(stats_df, "exact_match")

    inexact_mean_precision = compute_mean_precision(stats_df, "inexact_match")
    inexact_mean_recall = compute_mean_recall(stats_df, "inexact_match")

    return (
        exact_mean_precision,
        exact_mean_recall,
        inexact_mean_precision,
        inexact_mean_recall,
    )


def main():
    dir_path = Path(sys.argv[1])

    run_test = False
    if len(sys.argv) > 2:
        run_test = sys.argv[2] == "test"

    # load prefixes file
    prefixes = []
    with open(dir_path / "prefixes.txt", "r") as f:
        for line in f:
            prefixes.append(line.strip())

    overall_stats = []

    # run tool on each prefix
    for prefix in prefixes:
        # don't run test sets if not specified
        if "test" in prefix and not run_test:
            print(f"Skipping {prefix}")
            continue

        run_command(
            dir_path / prefix / f"{prefix}.repeat_seqs.tsv",
            dir_path / prefix / f"{prefix}.str_catalog.json",
            dir_path / prefix / f"{prefix}.strif_profile.tsv",
            visualize=True,
        )

        stat_df = create_stat_df(
            dir_path / prefix / f"{prefix}.truth.tsv",
            dir_path / prefix / f"{prefix}.strif_profile.tsv",
        )
        stat_df.to_csv(
            dir_path / prefix / f"{prefix}.compare.tsv", sep="\t", index=False
        )

        (
            exact_mean_prec,
            exact_mean_rec,
            inexact_mean_prec,
            inexact_mean_rec,
        ) = compute_metrics(stat_df)

        overall_stats.append(
            (
                prefix,
                exact_mean_prec,
                exact_mean_rec,
                inexact_mean_prec,
                inexact_mean_rec,
            )
        )

        print(f"Completed {prefix}:")
        print(
            f"  - Exact mean precision, recall: {exact_mean_prec:.2f}, {exact_mean_rec:.2f}"
        )
        print(
            f"  - Inexact mean precision, recall: {inexact_mean_prec:.2f}, {inexact_mean_rec:.2f}"
        )

    overall_stats_df = pd.DataFrame(
        overall_stats,
        columns=[
            "prefix",
            "exact_precision",
            "exact_recall",
            "inexact_precision",
            "inexact_recall",
        ],
    )
    overall_stats_df.to_csv(dir_path / f"overall_stats.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
