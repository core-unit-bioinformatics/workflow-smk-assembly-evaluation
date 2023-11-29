
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
        mem_mb = 16384
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
