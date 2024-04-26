#!/usr/bin/env python3

import argparse as argp
import collections as col
import enum
import hashlib as hl
import itertools as itt
import pathlib as pl
import sys

import pandas as pd

__PYTHON_VERSION = sys.version_info
if __PYTHON_VERSION.major < 3 and __PYTHON_VERSION.minor < 10:
    raise RuntimeError("Need at least Python v3.10 for execution")


Segment = col.namedtuple("Segment", ["start", "end", "label", "aln_context", "length", "block_id"])

Complement = col.namedtuple("Complement", ["start", "end", "label", "aln_context", "length", "block_id"])

class ComplementLabel(enum.Enum):
    UNALN = 1
    SCONTAIN = 2
    XCONTAIN = 3
    BRKOVL = 4
    ALNGAP = 5
    BRKSEQ = 6


def parse_command_line():

    parser = argp.ArgumentParser()

    parser.add_argument(
        "--alignments", "-a",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="alignment_file",
        help="Contig-to-reference alignments in normalized PAF format (TSV)."
    )

    parser.add_argument(
        "--target-label", "-tl",
        type=str,
        default="target",
        dest="target_label"
    )

    parser.add_argument(
        "--query-label", "-ql",
        type=str,
        default="query",
        dest="query_label"
    )

    parser.add_argument(
        "--delimiter", "--separator", "-d", "-s",
        type=str,
        choices=["@", "|", "_"],
        default="|",
        dest="delim"
    )

    parser.add_argument(
        "--out-target-all", "-ota",
        type=lambda x: pl.Path(x).resolve(strict=False),
        dest="out_target_all",
        default=None
    )

    parser.add_argument(
        "--out-target-complement", "-otc",
        type=lambda x: pl.Path(x).resolve(strict=False),
        dest="out_target_complement",
        default=None
    )

    parser.add_argument(
        "--out-query-all", "-oqa",
        type=lambda x: pl.Path(x).resolve(strict=False),
        dest="out_query_all",
        default=None
    )

    parser.add_argument(
        "--out-query-complement", "-oqc",
        type=lambda x: pl.Path(x).resolve(strict=False),
        default=None
    )

    parser.add_argument(
        "--dump-bed-like", "-bed",
        action="store_true",
        default=False,
        dest="dump_bed_like"
    )

    args = parser.parse_args()
    return args


def load_alignments(tsv_file):

    keep_columns = [
        "query_name", "query_length", "query_start", "query_end",
        "target_name", "target_length", "target_start", "target_end",
        "align_orient", "align_total", "mapq", "tp_align_type"
    ]

    aln = pd.read_csv(tsv_file, sep="\t", header=0, usecols=keep_columns)
    aln = aln.loc[aln["tp_align_type"] != 2, :].copy()

    return aln


def b_is_contained_in_a(segment_a, segment_b):
    """Function assumes that segments are processed
    in coordinate-ordered form,
    i.e. segment_a.start <= segment_b.start always holds.
    """
    return segment_b.start >= segment_a.start and segment_b.end <= segment_a.end


def a_overlaps_b(segment_a, segment_b):
    """Function assumes that is_contained(a,b)
    is checked first and segments are processed
    in coordinate-ordered sequence (see is_contained).
    Hence, overlap can only be forward-oriented.
    """
    return segment_b.start < segment_a.end


def build_segment_complements(sequence_length, segments, aln_context, delim):

    complements = []
    if segments[0].start > 0:
        # unaligned sequence at the beginning
        context_label = f"{aln_context}{delim}{segments[0].label}"
        complements.append(
            Complement(
                0, segments[0].start, ComplementLabel.UNALN.name,
                context_label, segments[0].start, f"START{delim}{segments[0].block_id}"
            )
        )

    for segment_a, segment_b in itt.pairwise(segments):
        nbh = f"{segment_a.block_id}{delim}{segment_b.block_id}"
        same_other_seq = segment_a.label == segment_b.label
        if b_is_contained_in_a(segment_a, segment_b):
            # by the way containment is checked,
            # segment_b must be smaller or equal
            if same_other_seq:
                label = ComplementLabel.SCONTAIN
                context_label = f"{aln_context}{delim}{segment_a.label}"
            else:
                label = ComplementLabel.XCONTAIN
                context_label = f"{aln_context}{delim}{segment_b.label}"
            complements.append(
                Complement(
                    segment_b.start, segment_b.end, label.name,
                    context_label, segment_b.end - segment_b.start,
                    nbh
                )
            )
            continue
        if a_overlaps_b(segment_a, segment_b):
            if same_other_seq:
                # ignore / alignment is collapsing for some reason
                continue
            # two non-identical seq. overlap, so there is
            # discontinuity in the query
            # by construction, this must be the overlap region
            ovl_start = segment_b.start
            ovl_end = segment_a.end
            assert ovl_start < ovl_end
            context_label = f"{aln_context}{delim}{segment_a.label}{delim}{segment_b.label}"
            complements.append(
                Complement(
                    ovl_start, ovl_end, ComplementLabel.BRKOVL.name,
                    context_label, ovl_end - ovl_start, nbh
                )
            )
            continue
        # subsequent alignments are interspersed
        if same_other_seq:
            # this could just be a break in the alignment or
            # indicate unresolved sequence / N-gap
            context_label = f"{aln_context}{delim}{segment_a.label}"
            gap_start, gap_end = segment_a.end, segment_b.start
            if gap_start == gap_end:
                gap_end += 1  # this is mainly to allow intersect ops.
            complements.append(
                Complement(
                    gap_start, gap_end, ComplementLabel.ALNGAP.name,
                    context_label, gap_end - gap_start, nbh
                )
            )
        else:
            # this should be a real break / gap in the query
            # because of two different (subsequent) query sequences
            context_label = f"{aln_context}{delim}{segment_a.label}{delim}{segment_b.label}"
            break_start, break_end = segment_a.end, segment_b.start
            if break_end == break_start:
                break_end += 1  # this is mainly to allow intersect ops.
            complements.append(
                Complement(
                    break_start, break_end, ComplementLabel.BRKSEQ.name,
                    context_label, break_end - break_start, nbh
                )
            )

    last_end = segments[-1].end
    last_seq = segments[-1].label
    # check if end of sequence is also unaligned
    if last_end < sequence_length:
        context_label = f"{aln_context}{delim}{last_seq}"
        complements.append(
            Complement(
                last_end, sequence_length, ComplementLabel.UNALN.name,
                context_label, sequence_length - last_end, f"{segments[-1].block_id}{delim}END"
            )
        )

    return complements


def determine_sequence_segments(alignments, target_label, query_label, delimiter):

    target_seq_sizes = dict()
    query_seq_sizes = dict()

    target_segments = col.defaultdict(list)
    query_segments = col.defaultdict(list)

    for row in alignments.itertuples(index=False):
        block_id = hl.md5("".join(map(str, list(row))).encode("utf-8")).hexdigest()

        target_seq_sizes[row.target_name] = row.target_length
        query_seq_sizes[row.query_name] = row.query_length

        # step 1: create segment in target coordinates -> use query_label
        segment_label = (
            f"{query_label}{delimiter}{row.query_name}{delimiter}"
            f"{row.query_start}-{row.query_end}{delimiter}{row.align_orient}"
            f"{delimiter}Q:{row.mapq}"
        )
        target_segments[row.target_name].append(
            Segment(row.target_start, row.target_end, row.query_name, segment_label, row.align_total, block_id)
        )

        # step 2: create segment in query coordinates -> use target_label
        segment_label = (
            f"{target_label}{delimiter}{row.target_name}{delimiter}"
            f"{row.target_start}-{row.target_end}{delimiter}{row.align_orient}"
            f"{delimiter}Q:{row.mapq}"
        )
        query_segments[row.query_name].append(
            Segment(row.query_start, row.query_end, row.target_name, segment_label, row.align_total, block_id)
        )

    return target_segments, target_seq_sizes, query_segments, query_seq_sizes


def build_dataframes(segments, seq_sizes, order_by, context_label, delim):

    all_regions = []
    complement_only = []
    for seq in seq_sizes.keys():
        seq_segments = sorted(segments[seq])
        df = pd.DataFrame.from_records(
            seq_segments, columns=["start", "end", "name", "aln_context", "length", "block_id"]
        )
        df["seq_name"] = seq
        df["seq_size"] = seq_sizes[seq]
        all_regions.append(df)
        seq_complements = build_segment_complements(seq_sizes[seq], seq_segments, context_label, delim)
        df = pd.DataFrame.from_records(
            seq_complements, columns=["start", "end", "name", "aln_context", "length", "block_id"]
        )
        df["seq_name"] = seq
        df["seq_size"] = seq_sizes[seq]
        all_regions.append(df)
        complement_only.append(df)

    all_regions = pd.concat(all_regions, axis=0, ignore_index=False)
    complement_only = pd.concat(complement_only, axis=0, ignore_index=False)

    if order_by == "lexorder":
        all_regions.sort_values(["seq_name", "start", "end"], inplace=True)
        complement_only.sort_values(["seq_name", "start", "end"], inplace=True)
    elif order_by == "size":
        all_regions.sort_values(
            ["seq_size", "seq_name", "start", "end"],
            ascending=[False, True, True, True], inplace=True
        )
        complement_only.sort_values(
            ["seq_size", "seq_name", "start", "end"],
            ascending=[False, True, True, True], inplace=True
        )
    else:
        raise ValueError(order_by)

    all_regions.reset_index(drop=True, inplace=True)
    all_regions = all_regions[["seq_name", "start", "end", "name", "aln_context", "length", "seq_size", "block_id"]]
    assert (all_regions["length"] > 0).all(), all_regions.loc[all_regions["length"] < 1, :]
    complement_only.reset_index(drop=True, inplace=True)
    complement_only = complement_only[["seq_name", "start", "end", "name", "aln_context", "length", "seq_size", "block_id"]]

    return all_regions, complement_only


def main():

    args = parse_command_line()

    alns = load_alignments(args.alignment_file)

    target_label = args.target_label.strip(args.delim)
    assert len(target_label) > 1
    query_label = args.query_label.strip(args.delim)
    assert len(query_label) > 1

    segment_info = determine_sequence_segments(alns, target_label, query_label, args.delim)
    target_segments, target_seq_sizes = segment_info[:2]
    query_segments, query_seq_sizes = segment_info[2:]

    target_regions, target_complements = build_dataframes(
        target_segments, target_seq_sizes, "lexorder",
        query_label, args.delim
    )

    if args.out_target_all is not None:
        args.out_target_all.parent.mkdir(exist_ok=True, parents=True)
        if args.dump_bed_like:
            target_regions.rename({"seq_name": "#seq_name"}, axis=1, inplace=True)
        target_regions.to_csv(args.out_target_all, sep="\t", header=True, index=False)

    if args.out_target_complement is not None:
        args.out_target_complement.parent.mkdir(exist_ok=True, parents=True)
        if args.dump_bed_like:
            target_complements.rename({"seq_name": "#seq_name"}, axis=1, inplace=True)
        target_complements.to_csv(args.out_target_complement, sep="\t", header=True, index=False)

    query_regions, query_complements = build_dataframes(
        query_segments, query_seq_sizes, "size",
        target_label, args.delim
    )

    if args.out_query_all is not None:
        args.out_query_all.parent.mkdir(exist_ok=True, parents=True)
        if args.dump_bed_like:
            query_regions.rename({"seq_name": "#seq_name"}, axis=1, inplace=True)
        query_regions.to_csv(args.out_query_all, sep="\t", header=True, index=False)

    if args.out_query_complement is not None:
        args.out_query_complement.parent.mkdir(exist_ok=True, parents=True)
        if args.dump_bed_like:
            query_complements.rename({"seq_name": "#seq_name"}, axis=1, inplace=True)
        query_complements.to_csv(args.out_query_complement, sep="\t", header=True, index=False)

    return 0


if __name__ == "__main__":
    main()
