#!/usr/bin/env python3

import argparse as argp
import contextlib as ctl
import pathlib as pl
import sys

import pandas as pd
import dnaio


def parse_command_line():

    parser = argp.ArgumentParser()

    parser.add_argument(
        "--input",
        "--fasta",
        "-i",
        "-f",
        dest="input_fasta",
        type=lambda x: pl.Path(x).resolve(strict=True),
        help="Path to FCS-processed FASTA file, with contigs flagged with source."
    )

    parser.add_argument(
        "--adapter-table",
        "-at",
        dest="adapter_table",
        type=lambda x: pl.Path(x).resolve(strict=True),
        help="Path to normalized adapter contamination report in TSV format."
    )

    parser.add_argument(
        "--contamination-table",
        "-ct",
        dest="contam_table",
        type=lambda x: pl.Path(x).resolve(strict=True),
        help="Path to normalized foreign contaminants report in TSV format."
    )

    parser.add_argument(
        "--out-pattern",
        "-op",
        dest="out_pattern",
        type=str,
        help="Split input FASTA by tags by adding the keyword SEQTAG in the output file pattern."
    )

    parser.add_argument(
        "--out-contaminants",
        "-oc",
        dest="out_contam",
        type=lambda x: pl.Path(x).resolve(),
        help="Path to FASTA output file for excluded sequences."
    )

    parser.add_argument(
        "--filter-tags",
        "-ft",
        dest="filter_tags",
        type=str,
        default=["hap1", "hap2", "unassigned", "disconnected", "mito", "rdna"],
        nargs="+",
        help="For this set of sequence tags, put contaminated sequences in the '--out-contaminants' file."
    )

    parser.add_argument(
        "--strip-tags",
        "-st",
        action="store_true",
        default=False,
        help="Remove sequence tags in output. Default: False"
    )

    parser.add_argument(
        "--report",
        "-r",
        action="store_true",
        default=False,
        dest="report",
        help="If true, write a brief summary report to stderr. Default: False"
    )

    args = parser.parse_args()
    return args


def read_screening_reports(args):

    adaptor = pd.read_csv(args.adapter_table, sep="\t", comment="#")
    adaptor_pass = set(adaptor.loc[adaptor["action"] == "PASS", "name"].values)

    contam = pd.read_csv(args.contam_table, sep="\t", comment="#")
    contam_pass = set(contam.loc[contam["action"] == "PASS", "name"].values)

    pass_sequences = contam_pass.intersection(adaptor_pass)

    return pass_sequences, adaptor, contam


def get_splitfile(out_pattern, seqtag):

    splitfile_name = out_pattern.replace("SEQTAG", seqtag)
    splitfile = pl.Path(splitfile_name).resolve()
    splitfile.parent.mkdir(exist_ok=True, parents=True)
    return splitfile


def get_normalized_entity_name(report, name, name_column):

    entity_name = report.loc[report["name"] == name, name_column].values[0]
    entity_name = entity_name.strip('"').replace(" ", ".")
    assert " " not in entity_name
    return entity_name


def trim_contaminated_sequence(report, name, lookup_name, sequence, name_column):

    trim_start, trim_end = report.loc[report["name"] == lookup_name, ["action_start", "action_end"]]
    seq_length = report.loc[report["name"] == lookup_name, "seq_length"].values[0]
    assert trim_start == 0 or trim_end == seq_length
    add_name = get_normalized_entity_name(report, name_column)

    if trim_start == 0:
        discard_seq = sequence[:trim_end]
        trimmed_seq = sequence[trim_end:]
    elif trim_end == seq_length:
        discard_seq = sequence[trim_start:]
        trimmed_seq = sequence[:trim_start]

    discard_header = f"{name}|TRIM|{trim_start}-{trim_end}|{add_name}"

    return trimmed_seq, discard_header, discard_seq


def process_contaminated_sequence(name, seqtag, sequence, adaptor_report, contam_report, filter_tags):

    lookup_name = name
    try:
        contam_action = contam_report.loc[contam_report["name"] == lookup_name, "action"].values[0]
    except IndexError:
        lookup_name = f"{name}.{seqtag}"
        contam_action = contam_report.loc[contam_report["name"] == lookup_name, "action"].values[0]
    adaptor_action = adaptor_report.loc[adaptor_report["name"] == lookup_name, "action"].values[0]

    # EXCLUDE flagged sequences are simply removed
    # from the main sequence output files specified
    # by the list of respective tags
    discard = seqtag in filter_tags
    if contam_action != "PASS":
        # contamination report has priority over adaptor report
        if contam_action == "EXCLUDE":
            trimmed_seq = None
            trimmed_header = None
            discard_seq = sequence
            discard_reason = get_normalized_entity_name(contam_report, lookup_name, "tax_division")
            discard_header = f"{name}|EXCLUDE|contam|{discard_reason}"
        elif contam_action == "TRIM":
            trimmed_seq, discard_seq, discard_header = trim_contaminated_sequence(
                contam_report, name, lookup_name, sequence, "tax_division"
            )
            trimmed_header = name
            discard = False
        else:
            raise ValueError(f"Cannot process foreign contamination action: {contam_action} / {name} / {seqtag}")

    elif adaptor_action != "PASS":
        # contamination has a PASS, which implies
        # some residual adaptor sequence somewhere
        if adaptor_action == "EXCLUDE":
            # quite unlikely for adaptor contamination
            trimmed_seq = None
            trimmed_header = None
            discard_seq = sequence
            discard_reason = get_normalized_entity_name(contam_report, lookup_name, "adaptor_name")
            discard_header = f"{name}|EXCLUDE|contam|{discard_reason}"
        elif adaptor_action == "TRIM":
            trimmed_seq, discard_seq, discard_header = trim_contaminated_sequence(
                contam_report, name, lookup_name, sequence, "adaptor_name"
            )
            trimmed_header = name
            discard = False
        else:
            raise ValueError(f"Cannot process adaptor contamination action: {contam_action} / {name} / {seqtag}")
    else:
        raise ValueError(f"Cannot process contaminated sequence record: {name} / {seqtag}")

    return discard, trimmed_header, trimmed_seq, discard_header, discard_seq


def main():

    args = parse_command_line()

    if args.out_pattern in ["stdout", "-", "", "/dev/stdout", "out"]:
        onefile = sys.stdout.buffer
        splitfiles = None
    elif "SEQTAG" not in args.out_pattern:
        onefile = pl.Path(args.out_pattern)
        onefile.parent.mkdir(exist_ok=True, parents=True)
        splitfiles = None
    else:
        onefile = None
        splitfiles = dict()

    args.out_contam.parent.mkdir(exist_ok=True, parents=True)
    contam_out = args.out_contam

    pass_sequences, adaptor_report, contam_report = read_screening_reports(args)

    count_records_in = 0
    count_records_out = 0
    count_contam_out = 0
    count_trim_op = 0

    with ctl.ExitStack() as exs:
        if onefile is not None:
            onefile = exs.enter_context(dnaio.open(onefile, mode="w", fileformat="fasta"))

        contam_out = exs.enter_context(dnaio.open(contam_out, mode="w", fileformat="fasta"))

        input_fasta = exs.enter_context(dnaio.open(args.input_fasta, mode="r"))

        for record in input_fasta:
            count_records_in += 1
            stripped_name, seqtag = record.name.rsplit(".", 1)
            if args.strip_tags:
                out_name = stripped_name
            else:
                out_name = record.name
            if record.name in pass_sequences:
                if onefile is not None:
                    write_out = onefile
                else:
                    try:
                        write_out = splitfiles[seqtag]
                    except KeyError:
                        splitfile = get_splitfile(args.out_pattern, seqtag)
                        write_out = exs.enter_context(dnaio.open(splitfile, mode="w", fileformat="fasta"))
                        splitfiles[seqtag] = write_out
                count_records_out += 1
                write_out.write(out_name, record.sequence)
            else:
                discard, trim_header, trim_seq, discard_header, discard_seq = process_contaminated_sequence(
                    out_name, seqtag, record.sequence,
                    adaptor_report, contam_report,
                    args.filter_tags
                )
                if discard:
                    # write discard_seq entry to contaminants file
                    count_contam_out += 1
                    contam_out.write(discard_header, discard_seq)
                else:
                    if onefile is not None:
                        write_out = onefile
                    else:
                        try:
                            write_out = splitfiles[seqtag]
                        except KeyError:
                            splitfile = get_splitfile(args.out_pattern, seqtag)
                            write_out = exs.enter_context(dnaio.open(splitfile, mode="w", fileformat="fasta"))
                            splitfiles[seqtag] = write_out

                    if trim_header is None:
                        # sequence was not trimmed -> exclude,
                        # but is kept based on sequence tag,
                        # so write original entry to normal output
                        count_records_out += 1
                        write_out.write(out_name, record.sequence)
                    else:
                        assert discard_header is not None
                        # sequence was trimmed, so write trimmed
                        # sequence file to output, and discard seq.
                        # to contaminants
                        count_records_out += 1
                        write_out.write(trim_header, trim_seq)
                        count_contam_out += 1
                        contam_out.write(discard_header, discard_seq)
                        count_trim_op += 1

    assert count_records_out + count_contam_out >= count_records_in

    if args.report:
        summary = (
            f"\n\n=== filter_contaminants report ===",
            f"\nRecords read from input: {count_records_in}"
            f"\nRecords written to output: {count_records_out}"
            f"\nRecords separated as contaminated: {count_contam_out}"
            f"\nRecords with trimmed sequences: {count_trim_op}\n\n"
        )
        sys.stderr.write(summary)

    return 0


if __name__ == "__main__":
    sys.exit(main())
