
def get_repeatmasker_run_memory_mb(input_size_mb, compressed=False):

    threshold_tiny = 1
    threshold_small = 100
    threshold_normal = 1000

    if compressed:
        compression_scaling = 1
    else:
        compression_scaling = 4  # ~gzip compressed FASTA vs uncompressed

    if input_size_mb < threshold_tiny * compression_scaling:
        mem_mb = 4096
    elif input_size_mb < threshold_small * compression_scaling:
        mem_mb = 49152
    elif input_size_mb < threshold_normal * compression_scaling:
        mem_mb = 98304
    else:
        mem_mb = 163840
    return mem_mb


def get_repeatmasker_run_time_hrs(input_size_mb, compressed=False):

    threshold_tiny = 1
    threshold_small = 100
    threshold_normal = 1000

    if compressed:
        compression_scaling = 1
    else:
        compression_scaling = 4  # ~gzip compressed FASTA vs uncompressed

    if input_size_mb < threshold_tiny * compression_scaling:
        time_hrs = 0
    elif input_size_mb < threshold_small * compression_scaling:
        time_hrs = 1
    elif input_size_mb < threshold_normal * compression_scaling:
        time_hrs = 71
    else:
        time_hrs = 84
    return time_hrs


def hmmer_scaling(resource, motif_name):

    assert resource in ["cpu", "mem", "time"]
    try:
        scaling_factor = int(HMMER_MOTIF_SEARCH[motif_name][f"scale_{resource}"])
    except KeyError:
        scaling_factor = 1
    return scaling_factor


def hmmer_threshold_value(threshold, motif_name):

    # defaults as per HMMER cli interface
    _DEFAULT_HMMER_EVALUE = "10"
    _DEFAULT_HMMER_SCORE = 0

    assert threshold in ["evalue_t", "evalue", "score", "score_t"]

    t_value = None
    if threshold in ["evalue", "evalue_t"]:
        # NB: only the E-value threshold is used in the HMMER call
        t_value = HMMER_MOTIF_SEARCH[motif_name].get("evalue_t", _DEFAULT_HMMER_EVALUE)

    if threshold in ["score", "score_t"]:
        t_value = HMMER_MOTIF_SEARCH[motif_name].get("score_t", _DEFAULT_HMMER_EVALUE)

    assert t_value is not None

    return t_value
