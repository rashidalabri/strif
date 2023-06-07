import os
from pathlib import Path
import sys
import pandas as pd
import random
import json


ALPHABET = ["A", "C", "G", "T"]
SEED = 42

DEFAULT_RANGES = {
    "motif": (2, 6),
    "seq": (10, 100),
    "intrpt": (0, 6),
}

N_SMALL = 10
N_LARGE = 1000
N_XLARGE = 10000


def random_seq(length, alphabet, allow_homopolymer=True):
    """
    Generate a random sequence of length n with given alphabet.
    Ensures that resulting sequence is not a homopolymer.
    """
    if length >= 2 and not allow_homopolymer:
        motif = random.sample(alphabet, 2) + random.choices(alphabet, k=length - 2)
    else:
        motif = random.choices(alphabet, k=length)
    return "".join(motif)


def rotate_seq(seq, n):
    """
    Rotate sequence by `n` positions.
    """
    return seq[n:] + seq[:n]


def repeat_seq(motif, n, rotate=True):
    """
    Generate a repeat sequence of length `n`.
    """
    n_repeat = n // len(motif) + 1
    seq = motif * n_repeat
    if rotate:
        rotate_n = random.randint(0, len(motif) - 1)
        seq = rotate_seq(seq, rotate_n)
    return seq[:n]


def simulate_repeat_seq(
    motif_len_range=DEFAULT_RANGES["motif"],
    seq_len_range=DEFAULT_RANGES["seq"],
    intrpt_len_range=DEFAULT_RANGES["intrpt"],
    insert=None,
    intersect_alpha=None,
    rotate=True,
):
    motif_len = random.randint(*motif_len_range)
    seq_len = random.randint(*seq_len_range)
    intrpt_len = random.randint(*intrpt_len_range)

    # if true,
    if intersect_alpha is None or intersect_alpha:
        motif_alpha = ALPHABET
        intrpt_alpha = ALPHABET
    else:
        # the motif and interruption sequences are generated from disjoint alphabets
        motif_alpha_len = random.randint(2, len(ALPHABET) - 1)
        motif_alpha = random.sample(ALPHABET, motif_alpha_len)
        intrpt_alpha = [x for x in ALPHABET if x not in motif_alpha]

    # generate repeat sequence
    motif = random_seq(motif_len, motif_alpha, allow_homopolymer=False)
    seq = repeat_seq(motif, seq_len, rotate=rotate)

    if intersect_alpha is not None and intersect_alpha:
        intrpt_alpha = list(set(motif))

    # generate interruption sequence and position
    intrpt = random_seq(intrpt_len, intrpt_alpha)
    intrpt_pos = random.randint(1, len(seq) - len(intrpt) - 1)

    # if insert is not specified, randomly choose whether to insert or substitute
    if insert is None:
        insert = bool(random.getrandbits(1))

    if insert:
        # insert the interruption sequence into the repeat sequence
        intrpt_seq = seq[:intrpt_pos] + intrpt + seq[intrpt_pos:]
        intrpt_seq = intrpt_seq[:seq_len]
    else:
        # substitute the interruption sequence for a portion of the repeat sequence
        intrpt_seq = seq[:intrpt_pos] + intrpt + seq[len(intrpt) + intrpt_pos :]

    if not insert and intrpt == seq[intrpt_pos : intrpt_pos + len(intrpt)]:
        # if the interruption sequence is the same as sequence
        # it is substituting, set the interruption sequence to be empty
        # since there will not be an interruption technically
        intrpt = ""

    return motif, intrpt, intrpt_seq


def generate_files(repeat_seqs, dir_path, prefix):
    motifs = [x[0] for x in repeat_seqs]
    intrpts = [x[1] for x in repeat_seqs]
    seqs = [x[2] for x in repeat_seqs]

    # create output directory
    os.makedirs(dir_path / prefix, exist_ok=True)

    # create file paths
    truth_path = dir_path / prefix / f"{prefix}.truth.tsv"
    repeat_seqs_path = dir_path / prefix / f"{prefix}.repeat_seqs.tsv"
    str_catalog_path = dir_path / prefix / f"{prefix}.str_catalog.json"

    # create truth dataframe
    df = pd.DataFrame(
        {
            "locus_id": range(len(repeat_seqs)),
            "motif": motifs,
            "interruption": intrpts,
            "seq": seqs,
        }
    )

    # write truth file
    df.to_csv(truth_path, index=False, sep="\t")

    # write repeat_seqs file for tool input
    df[["locus_id", "seq"]].to_csv(repeat_seqs_path, index=False, sep="\t", header=None)

    # create str_catalog file for tool input
    str_catalog = []
    for _, row in df.iterrows():
        motif = row["motif"]
        str_catalog.append(
            {"LocusId": str(row["locus_id"]), "LocusStructure": f"({motif})*"}
        )

    # write str_catalog file
    with open(str_catalog_path, "w") as f:
        json.dump(str_catalog, f, indent=4)


def generate_dataset(n, prefix, prefixes, dir_path, **kwargs):
    repeat_seqs = [simulate_repeat_seq(**kwargs) for _ in range(n)]
    generate_files(repeat_seqs, dir_path, prefix)
    prefixes.append(prefix)
    print(f"Generated {prefix} dataset")


def main():
    random.seed(SEED)
    dir_path = Path(sys.argv[1])

    prefixes = []

    # simple dataset
    generate_dataset(
        N_SMALL,
        "simple",
        prefixes,
        dir_path,
        motif_len_range=(2, 2),
        seq_len_range=(10, 10),
        intrpt_len_range=(1, 1),
        intersect_alpha=False,
        insert=True,
        rotate=False,
    )

    # no interruption dataset
    generate_dataset(
        N_SMALL,
        "no_interruption",
        prefixes,
        dir_path,
        intrpt_len_range=(0, 0),
    )

    # disjoint alphabet validation set (with variable interruption length)
    for i in range(1, 7):
        generate_dataset(
            N_LARGE,
            f"disjoint_{i}",
            prefixes,
            dir_path,
            intrpt_len_range=(i, i),
            intersect_alpha=False,
        )

    # intersecting alphabet validation set (with variable interruption length)
    for i in range(1, 7):
        generate_dataset(
            N_LARGE,
            f"intersect_{i}",
            prefixes,
            dir_path,
            intrpt_len_range=(i, i),
            intersect_alpha=True,
        )

    # insertion validation set
    for i in range(1, 7):
        generate_dataset(
            N_LARGE,
            f"insert_{i}",
            prefixes,
            dir_path,
            intrpt_len_range=(i, i),
            insert=True,
        )

    # substitution validation set
    for i in range(1, 7):
        generate_dataset(
            N_LARGE,
            f"substitute_{i}",
            prefixes,
            dir_path,
            intrpt_len_range=(i, i),
            insert=False,
        )

    # basic dataset
    for i in range(1, 7):
        generate_dataset(
            N_LARGE,
            f"basic_{i}",
            prefixes,
            dir_path,
            intrpt_len_range=(i, i),
        )

    # comprehensive set for training
    generate_dataset(
        N_LARGE,
        "comprehensive_train",
        prefixes,
        dir_path,
    )

    # comprehensive set for validation
    generate_dataset(
        N_LARGE,
        "comprehensive_valid",
        prefixes,
        dir_path,
    )

    # comprehensive set for testing
    generate_dataset(
        N_LARGE,
        "comprehensive_test",
        prefixes,
        dir_path,
    )

    # write prefixes file
    with open(dir_path / "prefixes.txt", "w") as f:
        f.write("\n".join(prefixes))


if __name__ == "__main__":
    main()
