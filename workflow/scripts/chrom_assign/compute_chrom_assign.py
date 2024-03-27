#!/usr/bin/env python3

import argparse as argp
import io
import pathlib as pl
import sys

import pandas as pd
import xopen


def parse_command_line():

    parser = argp.ArgumentParser()
    parser.add_argument(
        "--input",
        "-i",
        "-in-tsv",
        dest="input",
        type=lambda x: pl.Path(x).resolve(strict=True),
        help="Path to normalized PAF file as TSV.",
    )

    parser.add_argument(
        "--out-query",
        "-oq",
        dest="output_query",
        type=lambda x: pl.Path(x).resolve(strict=False),
        help="Path query-centric TSV output table."
    )

    parser.add_argument(
        "--out-target",
        "-ot",
        dest="output_target",
        type=lambda x: pl.Path(x).resolve(strict=False),
        help="Path target-centric TSV output table."
    )

    args = parser.parse_args()

    return args


def mle_target_seq_assignment(alignments):
    pass


def summarize_assignments(alignments):

    group_by = ["query_name", "target_name"]

    agg_stats = alignments.groupby(group_by).agg(
        query_length=("query_length", max),
        target_length=("target_length", max),
        align_total=("align_total", sum),
        align_matching=("align_matching", sum),
        align_longest=("align_matching", max),
        align_num=("align_matching", "nunique"),
        weights_mean=("weights_combined", "mean"),
        weights_median=("weights_combined", "median"),
        weights_max=("weights_combined", "max")
    )

    agg_stats["match_cov_query"] = (agg_stats["align_matching"] / agg_stats["query_length"]).round(5)
    agg_stats["match_cov_target"] = (agg_stats["align_matching"] / agg_stats["target_length"]).round(5)
    agg_stats["query_target_ratio"] = (agg_stats["query_length"] / agg_stats["target_length"]).round(5)
    agg_stats.reset_index(drop=False, inplace=True)

    return agg_stats


def compute_assignment_weights(alignments):


    alignments["weight_length"] = alignments["align_total"] / alignments["query_length"]
    alignments["weight_mapq"] = alignments["mapq"] / alignments["mapq"].max()
    alignments["weight_matched"] = alignments["align_matching"] / alignments["align_total"]
    alignments["weights_combined"] = (
        alignments["weight_length"]
        * alignments["weight_mapq"]
        * alignments["weight_matched"]
    )
    return alignments


def write_output(assignments, column_order, outfile, outheader):

    outfile.parent.mkdir(exist_ok=True, parents=True)

    with xopen.xopen(outfile, "w") as table:
        _ = table.write(outheader.getvalue())
        assignments[column_order].to_csv(
            table, sep="\t",
            header=True, index=False
        )
    return


def load_table_into_buffer(table_file):
    """This function exists for the annoying condition
    that someone used the hash '#' inside of a FASTA
    header as a delimiter, which cannot be handled
    by Pandas' comment policy in read_csv:
    > must be single character, stops parsing anywhere in the line

    Hence, this function loads the entire table into a buffer
    while ignoring only lines that start with a hash #.

    Args:
        table_file (pathlib.Path): File path to alignment table
    """

    table_buffer = io.StringIO()
    with xopen.xopen(table_file, "r") as table:
        for line in table:
            if line.startswith("#"):
                continue
            if not line.strip():
                continue
            table_buffer.write(line)
    # important to reset the buffer to position 0,
    # otherwise pandas.read_csv() will result in
    # empty dataframe
    table_buffer.seek(0)
    return table_buffer


def main():

    args = parse_command_line()

    load_columns = [
        "query_name", "query_length",
        "target_name", "target_length",
        "mapq",
        "align_total", "align_matching",
        "tp_align_type"
    ]

    # 2024-03-27 the following to accommodate inputs that contain
    # identifiers that use # as separators - yes, they exist ...
    try:
        alignments = pd.read_csv(
            args.input, sep="\t", header=0,
            usecols=load_columns, comment="#"
        )
    except ValueError:
        err_msg = (
            f"\nError parsing file: {args.input}\n"
            "Standard parsing with pandas.read_csv(... comment='#' ...) failed.\n"
            "Assuming that '#' was used as part of an identifier in the file.\n"
            "Loading table into buffer and restart...\n\n"
        )
        sys.stderr.write(err_msg)
        table_buffer = load_table_into_buffer(args.input)
        alignments = pd.read_csv(
            table_buffer, sep="\t", header=0,
            usecols=load_columns
        )

    out_info_header = io.StringIO()
    out_info_header.write(
        f"# total alignments: {alignments.shape[0]}\n"
        f"# total query sequences: {alignments['query_name'].nunique()}\n"
        f"# total target sequences: {alignments['target_name'].nunique()}\n"
    )
    # drop secondary or unaligned records if present
    alignments = alignments.loc[alignments["tp_align_type"].abs() == 1, :].copy()
    out_info_header.write(
        f"# non-secondary alignments: {alignments.shape[0]}\n"
    )
    alignments = alignments.loc[alignments["mapq"] > 0, :].copy()
    out_info_header.write(
        f"# non-secondary and MAPQ > 0: {alignments.shape[0]}\n"
        f"# final: total query sequences: {alignments['query_name'].nunique()}\n"
        f"# final: total target sequences: {alignments['target_name'].nunique()}\n"
    )
    alignments = compute_assignment_weights(alignments)
    assignments = summarize_assignments(alignments)

    generic_columns = [
        "query_target_ratio",
        "query_length",
        "match_cov_query",
        "target_length",
        "match_cov_target",
        "align_num",
        "align_total",
        "align_matching",
        "weights_mean",
        "weights_median",
        "weights_max"
    ]

    sort_order_query = ["query_name", "target_name"] + generic_columns
    assignments.sort_values(
        ["query_name", "match_cov_query"],
        ascending=[True, False],
        inplace=True
    )
    write_output(
        assignments, sort_order_query,
        args.output_query, out_info_header
    )

    sort_order_target = ["target_name", "query_name"] + generic_columns
    assignments.sort_values(
        ["target_name", "match_cov_target"],
        ascending=[True, False],
        inplace=True
    )
    write_output(
        assignments, sort_order_target,
        args.output_target, out_info_header
    )

    return 0



if __name__ == "__main__":
    sys.exit(main())
