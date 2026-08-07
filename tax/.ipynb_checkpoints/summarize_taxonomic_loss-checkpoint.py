#!/usr/bin/env python3

"""
Compare taxonomic signal loss caused by assembly.

trim:
    original trimmed metagenome sourmash gather

inflate:
    assembly-recovered k-mers with original read abundances

Outputs:
    1. per-sample summary
    2. per-sample-family retention table
"""

from pathlib import Path
import pandas as pd
import numpy as np


# Titus outputs
TRIM_DIR = Path(
    "/home/ctbrown/scratch3/2026-assembly-loss-tax/outputs/trim"
)

INFLATE_DIR = Path(
    "/home/ctbrown/scratch3/2026-assembly-loss-tax/outputs/inflate"
)


# local outputs
OUT_DIR = Path(
    "/home/zyzhao/2025-zyzhao-assemloss/tax/outputs"
)

OUT_DIR.mkdir(exist_ok=True, parents=True)


def load_family_counts(filename):
    """
    Collapse sourmash gather results to family level.
    """

    df = pd.read_csv(filename)

    # GTDB lineage format:
    # d__;p__;c__;o__;f__;g__;s__
    df["family"] = (
        df["lineage"]
        .fillna("")
        .str.split(";")
        .str[4]
    )

    df = df[
        (df["family"].notna()) &
        (df["family"] != "")
    ]

    family_counts = (
        df.groupby("family")
        ["n_unique_weighted_found"]
        .sum()
    )

    return family_counts


sample_summary = []
family_summary = []


trim_files = sorted(
    TRIM_DIR.glob("*.gather.with-lineages.csv")
)


print(f"Found {len(trim_files)} trim files")


for trim_file in trim_files:

    sample = trim_file.name.replace(
        ".gather.with-lineages.csv",
        ""
    )

    inflate_file = (
        INFLATE_DIR /
        trim_file.name
    )

    if not inflate_file.exists():
        print(
            f"Missing inflate file for {sample}"
        )
        continue


    print(sample)


    trim = load_family_counts(trim_file)
    inflate = load_family_counts(inflate_file)


    combined = pd.concat(
        [
            trim.rename("trim"),
            inflate.rename("inflate")
        ],
        axis=1
    ).fillna(0)


    combined["retention"] = np.where(
        combined["trim"] > 0,
        combined["inflate"] /
        combined["trim"],
        np.nan
    )


    trim_total = combined["trim"].sum()
    inflate_total = combined["inflate"].sum()


    combined["trim_relative"] = (
        combined["trim"] /
        trim_total
    )

    combined["inflate_relative"] = (
        combined["inflate"] /
        inflate_total
    )


    combined["relative_shift"] = (
        combined["inflate_relative"]
        -
        combined["trim_relative"]
    )


    combined = (
        combined
        .reset_index()
        .rename(columns={"index": "family"})
    )


    combined.insert(
        0,
        "sample",
        sample
    )


    family_summary.append(
        combined
    )


    sample_summary.append(
        {
            "sample": sample,

            "trim_taxonomic_signal":
                trim_total,

            "inflate_taxonomic_signal":
                inflate_total,

            "overall_retention":
                inflate_total / trim_total,

            "trim_family_count":
                int(
                    (combined["trim"] > 0)
                    .sum()
                ),

            "inflate_family_count":
                int(
                    (combined["inflate"] > 0)
                    .sum()
                ),

            "families_lost":
                int(
                    (
                        (combined["trim"] > 0)
                        &
                        (combined["inflate"] == 0)
                    )
                    .sum()
                )
        }
    )


family_df = pd.concat(
    family_summary,
    ignore_index=True
)


sample_df = pd.DataFrame(
    sample_summary
)


family_df.to_csv(
    OUT_DIR /
    "family_retention_by_sample.csv",
    index=False
)


sample_df.to_csv(
    OUT_DIR /
    "assembly_taxonomic_retention_by_sample.csv",
    index=False
)


print()
print("Finished")
print(
    sample_df[
        [
            "sample",
            "overall_retention",
            "families_lost"
        ]
    ]
)
