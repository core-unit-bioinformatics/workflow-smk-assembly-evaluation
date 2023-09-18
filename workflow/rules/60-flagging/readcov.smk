
rule prepare_window_read_coverage_histogram:
    input:
        hdf = rules.transform_mosdepth_window_read_coverage.output.hdf
    output:
        tsv = DIR_PROC.joinpath(
            "60-flagging", "readcov", "window_histogram",
            "{sample}.read-cov-hist.tsv"
        )
    resources:
        mem_mb=lambda wildcards, attempt: 2048 * attempt
    params:
        low_bin=-1,
        high_bin=301,
        step_size=1,
        mapq=0,
        signal="pct_median_cov"
    run:
        import numpy as np
        import pandas as pd
        import pathlib as pl
        import sys

        bins = np.arange(params.low_bin, params.high_bin, step=params.step_size).round(0)
        bins = np.hstack((bins, np.array([params.high_bin, sys.maxsize])))

        bin_labels = []
        for edge_idx, bin_edge in enumerate(bins[1:], start=1):
            if edge_idx == bins.size - 1:
                label = "LRG"
            else:
                label = f"{bin_edge:03}"
            bin_labels.append(label)

        # TODO move to settings module
        cov_groups = config["groups_window_read_coverage"]

        group_counts = dict()
        with pd.HDFStore(input.hdf, "r") as hdf:
            for key in hdf.keys():
                use_group = None
                for group_name in cov_groups:
                    if key.startswith(f"/{group_name}"):
                        use_group = group_name
                        break
                if use_group is None:
                    continue

                data = hdf[key]
                for column in data.columns:
                    read_type, aln_type, mapq, stat = column
                    if stat != params.signal:
                        continue
                    if mapq != params.mapq:
                        continue

                    hist = pd.cut(
                        data.loc[:, column],
                        bins=bins,
                        include_lowest=False, right=True,
                        labels=bin_labels
                    ).value_counts()
                    store_key = (wildcards.sample, read_type, aln_type, mapq, use_group)
                    if store_key not in group_counts:
                        tmp = pd.Series(
                            np.zeros(bins.size-1, dtype=int),
                            index=bin_labels
                        )
                        group_counts[store_key] = tmp
                    group_counts[store_key] += hist

        sample_hist = pd.DataFrame.from_dict(group_counts, orient="index")
        sample_hist.index.rename(["sample", "read_type", "aln_type", "mapq", "asm_unit"], inplace=True)
        sample_hist.to_csv(output.tsv, sep="\t", header=True, index=True)
    # END OF RUN BLOCK


rule run_all_window_read_coverage_histograms:
    input:
        tsv = expand(
            rules.prepare_window_read_coverage_histogram.output.tsv,
            sample=SAMPLES
        )

