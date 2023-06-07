import errno
import os
import sys
from itertools import product
from pathlib import Path
import time

import pandas as pd
from metrics import compute_metrics, create_stat_df, run_command

CORES = 10

MATCH_SCORE_RANGE = list(range(1, 2))
MISMATCH_PENALTY_RANGE = list(range(1, 11))
GAP_OPEN_PENALTY_RANGE = list(range(1, 11))
GAP_EXTEND_PENALTY_RANGE = list(range(1, 11))


def start_process(id, params, valid_dir_path, prefix, tmp_dir_path, results_file):
    pid = os.fork()
    if pid == 0:
        print(f"[Process {id}] Starting...", flush=True)

        truth_path = valid_dir_path / prefix / f"{prefix}.truth.tsv"
        repeat_seqs_path = valid_dir_path / prefix / f"{prefix}.repeat_seqs.tsv"
        str_catalog_path = valid_dir_path / prefix / f"{prefix}.str_catalog.json"

        for i, (m, x, g, e) in enumerate(params):
            profile_path = tmp_dir_path / prefix / f"{prefix}_proc{id}_{i}.tmp"

            if i % 100 == 0:
                percent_complete = int(i / len(params) * 100)
                print(f"[Process {id}] {percent_complete}%", flush=True)

            cmd = run_command(
                repeat_seqs_path, str_catalog_path, profile_path, m, x, g, e
            )
            time.sleep(0.01)

            try:
                stat_df = create_stat_df(truth_path, profile_path)
                metrics = "\t".join([str(m) for m in compute_metrics(stat_df)])
                results_file.write(f"{m}\t{x}\t{g}\t{e}\t{metrics}\n")
                results_file.flush()
            except OSError as e:
                if e.errno != errno.ENOENT:
                    raise
                print(f"[Process {id}] Failed to run command: {cmd}", flush=True)

            # delete profile file
            try:
                os.remove(profile_path)
            except OSError as e:
                if e.errno != errno.ENOENT:
                    raise

        # when finished, exit the child process
        print(f"[Process {id}] Finished", flush=True)
        results_file.close()
        exit(0)
    else:
        return pid


def sort_results(unsorted_file_path, sorted_file_path):
    # sort the results using pandas
    print("Sorting results...", flush=True)
    results = pd.read_csv(unsorted_file_path, sep="\t")
    results["mean_recall"] = (results["exact_recall"] + results["inexact_recall"]) / 2
    results.sort_values("mean_recall", ascending=False, inplace=True)
    results.to_csv(sorted_file_path, sep="\t", index=False)

    results = results.head(10)

    # print average of each column
    print("----------------------------------------", flush=True)
    print("Avg. of top 10 parameters:", flush=True)
    print(results.mean(axis=0), flush=True)
    print("----------------------------------------", flush=True)

    # delete unsorted file
    os.remove(unsorted_file_path)


def perform_param_grid_search(params, valid_dir_path, prefix, tmp_dir_path, cores):
    print(f"Testing {len(params)} combinations using {cores} cores...")

    with open(
        valid_dir_path / prefix / f"{prefix}.param_search.unsorted.tsv", "w"
    ) as f:
        f.write(
            "match_score\tmismatch_penalty\tgap_open_penalty\tgap_extend_penalty\texact_precision\texact_recall\tinexact_precision\tinexact_recall\n"
        )
        f.flush()

        batch_size = len(params) // cores + 1

        for i in range(cores):
            start_idx = i * batch_size
            batch = params[start_idx : start_idx + batch_size]
            start_process(i + 1, batch, valid_dir_path, prefix, tmp_dir_path, f)

    for _ in range(cores):
        os.wait()

    sort_results(
        valid_dir_path / prefix / f"{prefix}.param_search.unsorted.tsv",
        valid_dir_path / prefix / f"{prefix}.param_search.tsv",
    )

    print("Done", flush=True)


def main():
    params = list(
        product(
            MATCH_SCORE_RANGE,
            MISMATCH_PENALTY_RANGE,
            GAP_OPEN_PENALTY_RANGE,
            GAP_EXTEND_PENALTY_RANGE,
        )
    )

    # remove params where gap_open_penalty < gap_extend_penalty
    params = [(m, x, g, e) for m, x, g, e in params if x > m and g > e]

    valid_dir_path = Path(sys.argv[1])
    tmp_dir_path = valid_dir_path
    prefix = sys.argv[2]
    perform_param_grid_search(params, valid_dir_path, prefix, tmp_dir_path, CORES)


if __name__ == "__main__":
    main()
