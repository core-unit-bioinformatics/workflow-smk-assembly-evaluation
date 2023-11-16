
MOSDEPTH_OUTPUT = []

if RUN_MOSDEPTH and RUN_MOSDEPTH_ASSM_REF_COV:

    MOSDEPTH_OUTPUT.extend(
        rules.mosdepth_assembly_reference_coverage_window.input.windowed
    )
    MOSDEPTH_OUTPUT.extend(
        rules.mosdepth_assembly_reference_coverage_window.input.merged
    )

if RUN_MOSDEPTH and RUN_MOSDEPTH_ASSM_READ_COV:

    MOSDEPTH_OUTPUT.extend(
        rules.run_all_mosdepth_assembly_read_coverage.input.check_files
    )
    MOSDEPTH_OUTPUT.extend(
        rules.run_all_mosdepth_assembly_read_coverage.input.transform_cov
    )

    MOSDEPTH_OUTPUT.extend(
        rules.run_all_mosdepth_coverage_stats.input.stats
    )
