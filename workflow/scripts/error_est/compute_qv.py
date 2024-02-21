#!/usr/bin/env python3

import argparse as argp
import collections as col
import math
import pathlib as pl

import pandas as pd
import pysam


def parse_command_line():

    parser = argp.ArgumentParser()

    parser.add_argument(
        "--input", "-i",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="input",
        help="Path to input files."
    )

    parser.add_argument(
        "--subtract", "-s",
        type=lambda x: pl.Path(x).resolve(strict=False),
        default=None,
        dest="subtract",
        help="BED file of regions to ignore when computing the QV estimate."
    )

    parser.add_argument(
        "--mode", "-m",
        type=str,
        choices=["vcf"],
        default="vcf",
        dest="mode",
        help="Select mode for QV estimation: [vcf]-based ... Default: vcf"
    )

    parser.add_argument(
        "--output", "-o",
        type=lambda x: pl.Path(x).resolve(strict=False),
        dest="output",
        help="Paht to output file."
    )

    args = parser.parse_args()
    return args


def compute_qv(num_errors, ref_size):
    """_summary_

    Args:
        num_errors (_type_): _description_
        ref_size (_type_): _description_

    Returns:
        _type_: _description_
    """
    p = num_errors / ref_size
    try:
        q = -10 * math.log10(p)
    except ValueError:
        return 99
    return int(round(q, 0))


def read_vcf_records(input_file):

    error_counts = col.Counter()

    with pysam.VariantFile(input_file) as vcf:
        contig_lens = dict(
            (contig.name, contig.length) for contig in
            vcf.header.contigs.itervalues()
        )

        for record in vcf:
            ref_allele = record.ref
            assert isinstance(ref_allele, str)
            alt_alleles = record.alts

            max_length_diff = 0
            for alt in alt_alleles:
                # diff in allele length is zero for SNV,
                # hence count those as 1 via max()
                length_diff = max(1, abs(len(ref_allele) - len(alt)))
                max_length_diff = max(length_diff, max_length_diff)

            error_counts[record.chrom] += max_length_diff

    return error_counts, contig_lens


def process_vcf_errors(error_counts, contig_lengths, subtract_lengths):

    out_records = []
    for contig, contig_len in contig_lengths.items():
        # NB: this is a Counter(), default return 0
        error_bp = error_counts[contig]
        adjusted_length = contig_len - subtract_lengths[contig]
        if error_bp < 1:
            qv_est = 99
        else:
            assert 0 < adjusted_length <= contig_len
            qv_est = compute_qv(error_bp, adjusted_length)
        out_records.append(
            (contig, contig_len, adjusted_length, error_bp, qv_est)
        )

    total_errors = sum(error_counts.values())
    total_length = sum(contig_lengths.values())
    total_subtract = sum(subtract_lengths.values())
    adjusted_total = total_length - total_subtract
    assert 0 < adjusted_total <= total_length
    wg_qv_est = compute_qv(total_errors, adjusted_total)
    out_records.append(
        ("genome", total_length, adjusted_total, total_errors, wg_qv_est)
    )
    out_records = sorted(out_records, key=lambda t: t[1], reverse=True)

    return out_records


def dump_by_contig_qv_estimates(out_records, out_file):

    out_file.parent.mkdir(exist_ok=True, parents=True)

    df = pd.DataFrame.from_records(
        out_records,
        columns=["seq_name", "seq_length", "adj_length", "num_errors", "qv"]
    )

    df.to_csv(out_file, sep="\t", header=True, index=False)

    return


def read_subtract_lengths(subtract_lengths):

    if subtract_lengths is not None:
        df = pd.read_csv(
            subtract_lengths, sep="\t", header=None,
            comment="#", usecols=[0,1,2]
        )
        df.columns = ["contig", "start", "end"]
        df["length"] = df["end"] - df["start"]
        subtract_lookup = dict(
            (k, v) for k, v in df.groupby("contig")["length"].sum().items()
        )
        # NB: col.Counter() returns 0 for non-ex keys
        subtract_lookup = col.Counter(subtract_lookup)
    else:
        subtract_lookup = col.Counter()

    return subtract_lookup


def run_vcf_based_qv_estimate(input_vcf, output_table, subtract_lengths):

    error_counts, contig_lengths = read_vcf_records(input_vcf)
    subtract_lengths = read_subtract_lengths(subtract_lengths)
    out_records = process_vcf_errors(error_counts, contig_lengths, subtract_lengths)
    dump_by_contig_qv_estimates(out_records, output_table)

    return


def main():

    args = parse_command_line()

    if args.mode == "vcf":
        run_vcf_based_qv_estimate(args.input, args.output, args.subtract)
    else:
        raise NotImplementedError(args.mode)

    return 0


if __name__ == "__main__":
    main()
