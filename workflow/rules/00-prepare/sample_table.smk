import pathlib
import hashlib
import collections
import re

import pandas

SAMPLES = None
SAMPLE_INFOS = None


def process_sample_sheet():

    SAMPLE_SHEET_FILE = pathlib.Path(config["samples"]).resolve(strict=True)
    SAMPLE_SHEET = pandas.read_csv(
        SAMPLE_SHEET_FILE,
        sep="\t",
        header=0
    )

    # step 1: each row is a sample,
    # just collect the input files
    sample_infos = collect_sample_infos(SAMPLE_SHEET)
    all_samples = sorted(sample_infos.keys())
    assert len(all_samples) == SAMPLE_SHEET.shape[0]

    global SAMPLES
    SAMPLES = all_samples

    global SAMPLE_INFOS
    SAMPLE_INFOS = sample_infos

    return


def collect_sample_infos(sample_sheet):
    """
    Iterate over rows and collect sample
    infos; currently supported:
    - assembly unit (asm_)
    - reads (reads_)
    - sex (sex: f/m/u/d/x)
    """

    sample_infos = collections.defaultdict(dict)

    sheet_columns = sample_sheet.columns
    sheet_columns = [
        c for c in sheet_columns if c not in ["sample", "sex"]
    ]

    for row in sample_sheet.itertuples():
        sample = row.sample
        try:
            sample_sex = normalize_sample_sex(row.sex)
        except AttributeError:
            sample_sex = "any"
        sample_infos[sample]["sex"] = sample_sex
        for column in sheet_columns:
            if columns.startswith("asm_") or column.startswith("reads_"):
                group_type, group_id = columns.split("_", 1)
                check_data_identifier(group_id)
                # three lists: file names, full path hashes,
                # path hash prefixes [path IDs]
                seq_file_infos = collect_sequence_input(getattr(row, column))
            else:
                sample_infos[sample][column] = getattr(row, column)
                continue

            if group_type == "asm":
                assert len(seq_file_infos[0]) == 1
                # support one file per assembly unit,
                # path hash is irrelevant
                seq_file_infos[-1] = [None]
                seq_file_infos[-2] = [None]

            for seqfile, phash, pathid in seq_file_infos:
                sample_infos[sample][(group_type, group_id, phash)] = seqfile
                all_key = group_type, group_id, "hash_all"
                if all_key not in sample_infos[sample]:
                    sample_infos[sample][all_key] = {phash}
                else:
                    sample_infos[sample][all_key].add(phash)
                sample_infos[sample][(group_type, group_id, pathid)] = seqfile
                all_key = group_type, group_id, "path_id_all"
                if all_key not in sample_infos[sample]:
                    sample_infos[sample][all_key] = {pathid}
                else:
                    sample_infos[sample][all_key].add(pathid)

    return sample_infos


def check_data_identifier(identifier):
    if re.match("[a-z0-9]+$", identifier) is None:
        raise ValueError(f"Invalid data identifier: {identifier}")
    return


def normalize_sample_sex(sample_sex):

    known_norm = {
        "f": "female",
        "m": "male",
        "u": "any",
        "d": "any",
        "x": "any"
    }

    norm_sex = known_norm.get(sample_sex.lower(), sample_sex.lower())

    if not norm_sex in ["male", "female", "any"]:
        raise ValueError(f"Cannot normalize sample sex: {sample_sex}")

    return norm_sex


def collect_input_files(sample_sheet):
    """
    The output of this function should
    be sufficient to run the workflow
    in single-sample mode
    """
    sample_input = collections.defaultdict(dict)
    trio_samples = set()
    sseq_samples = set()
    unphased_samples = set()

    for row in sample_sheet.itertuples():
        sample = row.sample
        hifi_input, hifi_hashes = collect_sequence_input(row.hifi)
        ont_input, ont_hashes = collect_sequence_input(row.ont)
        if row.target == "trio":
            assert hasattr(row, "hap1") and hasattr(row, "hap2"), (
                f"Trio-phasing a sample requires hap1/hap2 "
                "columns in sample sheet: {sample}"
            )
            trio_samples.add(sample)
            hap1_db = row.hap1
            assert hap1_db.endswith("meryl")
            assert pathlib.Path(hap1_db).resolve(strict=True).is_dir()
            sample_input[sample]["hap1"] = hap1_db

            hap2_db = row.hap2
            assert hap2_db.endswith("meryl")
            assert pathlib.Path(hap2_db).resolve(strict=True).is_dir()
            sample_input[sample]["hap2"] = hap2_db

            assert sample_input[sample]["hap1"] != sample_input[sample]["hap2"]

        elif row.target == "sseq":
            # NB: at the moment, phasing an assembly with
            # Strand-seq requires to complete the unphased assembly
            # first, hence the sample is always added to the
            # unphased_samples set as well. This may need to be
            # changed if both pipelines can be integrated.
            unphased_samples.add(sample)
            sseq_samples.add(sample)
            assert hasattr(row, "phasing_paths"), (
                "Strand-seq phasing of sample requires"
                f"phasing_paths column in sample sheet: {sample}"
            )
            phasing_paths_file = pathlib.Path(row.phasing_paths).resolve(strict=True)
            assert phasing_paths_file.is_file(), (
                f"Phasing paths entry is not a valid file: {phasing_paths_file}"
            )
            assert phasing_paths_file.suffix == ".gaf", (
                f"Phasing paths file must be GAF: {phasing_paths_file.name}"
            )
            sample_input[sample]["phasing_paths"] = phasing_paths_file

        elif row.target == "unphased":
            unphased_samples.add(sample)
        else:
            raise ValueError(row.target)
        sample_input[sample]["hifi"] = hifi_input
        sample_input[sample]["ont"] = ont_input

    return sample_input, trio_samples, sseq_samples, unphased_samples


def _read_input_files_from_fofn(fofn_path):
    """Read input file listing from
    file of file names
    TODO: candidate for inclusion in template
    """

    input_files = []
    with open(fofn_path, "r") as listing:
        for line in listing:
            if not line.strip():
                continue
            try:
                file_path = pathlib.Path(line.strip()).resolve(strict=True)
            except FileNotFoundError:
                try:
                    file_path = DATA_ROOT.joinpath(line.strip()).resolve(strict=True)
                except FileNotFoundError:
                    err_msg = "\nERROR\n"
                    err_msg += f"Cannot find file: {line.strip}\n"
                    err_msg += f"Data root is set to: {DATA_ROOT}\n"
                    sys.stderr.write(err_msg)
                    raise
            input_files.append(file_path)

    return sorted(input_files)


def subset_path(full_path):
    """This helper exists to reduce
    the absolute path to a file
    to just the file name and its
    parent.
    TODO: should be codified as part
    of the template utilities to improve
    infrastructure portability of active
    workflows
    """
    folder_name = full_path.parent.name
    file_name = full_path.name
    subset_path = f"{folder_name}/{file_name}"
    # if it so happens that the file resides
    # in a root-level location, strip off
    # leading slash
    return subset_path.strip("/")


def collect_sequence_input(path_spec):
    """
    Generic function to collect HiFi or
    ONT/Nanopore input (read) files
    """
    input_files = []
    input_hashes = []
    input_path_ids = []
    for sub_input in path_spec.split(","):
        input_path = pathlib.Path(sub_input).resolve()
        if input_path.is_file() and input_path.name.endswith(".fofn"):
            fofn_files = _read_input_files_from_fofn(input_path)
            fofn_hashes = [
                hashlib.sha256(
                    subset_path(fp).encode("utf-8")
                ).hexdigest() for fp in fofn_files
            ]
            input_files.extend(fofn_files)
            input_hashes.extend(fofn_hashes)
            input_path_ids.extend(
                [fh[:PATH_ID_LENGTH] for fh in fofn_hashes]
            )
        elif input_path.is_file():
            input_hash = hashlib.sha256(
                subset_path(input_path).encode("utf-8")
            ).hexdigest()
            input_files.append(input_path)
            input_hashes.append(input_hash)
            input_path_ids.appnd(input_hash[:PATH_ID_LENGTH])
        elif input_path.is_dir():
            collected_files = _collect_files(input_path)
            collected_hashes = [
                hashlib.sha256(
                    subset_path(fp).encode("utf-8")
                ).hexdigest() for fp in collected_files
            ]
            input_files.extend(collected_files)
            input_hashes.extend(collected_hashes)
            input_path_ids.extend(
                [ch[:PATH_ID_LENGTH] for ch in collected_hashes]
            )
        else:
            raise ValueError(f"Cannot handle input: {sub_input}")
    return input_files, input_hashes, input_path_ids


def _build_constraint(values):
    escaped_values = sorted(map(re.escape, map(str, values)))
    constraint = "(" + "|".join(escaped_values) + ")"
    return constraint


def _collect_files(folder):

    all_files = set()
    for pattern in config["input_file_ext"]:
        pattern_files = set(folder.glob(f"**/*.{pattern}"))
        all_files = all_files.union(pattern_files)
    all_files = [f for f in sorted(all_files) if f.is_file()]
    if len(all_files) < 1:
        raise ValueError(f"No input files found underneath {folder}")
    return all_files

process_sample_sheet()
