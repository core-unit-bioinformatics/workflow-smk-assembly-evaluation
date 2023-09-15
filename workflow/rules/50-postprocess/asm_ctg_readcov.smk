
rule mosdepth_assembly_read_coverage_window:
    """Note that this rule assumes a pre-filtered
    BAM as input file. Thus, the exlude flag option
    is explicitly set to 0 (unset).
    """
    input:
        bam = rules.merge_read_to_assembly_subset_alignments.output.bam,
        bai = rules.merge_read_to_assembly_subset_alignments.output.bai
    output:
        check = DIR_PROC.joinpath(
            "50-postprocess", "asm_ctg_readcov", "mosdepth",
            "{sample}.{read_type}.{aln_subset}.mq{mapq}.wd",
            "{sample}.{read_type}.{aln_subset}.mq{mapq}.ok"
        )
    threads: CPU_LOW
    conda: DIR_ENVS.joinpath("biotools", "mosdepth.yaml")
    resources:
        mem_mb = lambda wildcards, attempt: 8192 * attempt,
        time_hrs = lambda wildcards, attempt: max(0, attempt-1)
    params:
        window_size = MOSDEPTH_ASSM_READ_COV_WINDOW_SIZE,
        min_mapq = lambda wildcards: int(wildcards.mapq),
        out_prefix = lambda wildcards, output: pathlib.Path(output.check).with_suffix("")
    shell:
        "mosdepth --threads {threads} --by {params.window_size} "
            "--no-per-base --flag 0 --use-median "
            "--quantize 0:5:8:12:15:20:25:30:50:100: --fast-mode --mapq {params.min_mapq} "
            "{params.out_prefix} {input.bam}"
            " && "
        "touch {output.check}"


rule mosdepth_coverage_stats_summary:
    """ Generate a summary info
    file of the (mean) coverages
    along the assemble sequences
    to be used for normalization
    in later stages

    Makes only sense for primary
    alignments at MAPQ 0
    """
    input:
        check_file = expand(
            rules.mosdepth_assembly_read_coverage_window.output.check,
            aln_subset=["onlyPRI"],
            mapq=["00"],
            allow_missing=True
        )
    output:
        stats = DIR_PROC.joinpath(
            "50-postprocess", "asm_ctg_readcov", "mosdepth",
            "{sample}.{read_type}.cov-stats.tsv"
        )
    params:
        threshold_percentile = 99.5
    run:
        import pathlib as pl
        import pandas as pd
        import numpy as np
        import scipy.stats as stats
        _this = "50-postprocess::asm_ctg_readcov::summarize_mosdepth_coverage"

        summary_file = pl.Path(input.check_file).with_suffix(".mosdepth.summary.txt")
        if not summary_file.is_file():
            logerr(f"\nERROR in {_this} - file does not exist {summary_file}\n")
            raise FileNotFoundError(summary_file.name)
        df = pd.read_csv(summary_file, sep="\t", header=0)
        df = df.loc[~df["chrom"].str.endswith("region"), :].copy()
        reported_total_mean = df.loc[df["chrom"] == "total", "mean"].values[0]
        df = df.loc[~df["chrom"].str.contains("total"), :].copy()

        wt_avg = np.average(df["mean"].values, weights=df["length"].values)
        wt_var = np.average((df["mean"].values - wt_avg)**2, weights=df["length"].values)
        wt_stddev = np.sqrt(wt_var)

        threshold_score = stats.scoreatpercentile(
            df["mean"].values, params.threshold_percentile
        )
        assert np.isclose(reported_total_mean, wt_avg, atol=0.5), \
            f"Averages do not match: {reported_total_mean} vs {wt_avg}"

        fmt_pctile = str(params.threshold_percentile).replace(".", "-")
        df = pd.DataFrame(
            [[
                wildcards.sample, wildcards.read_type, reported_total_mean,
                wt_avg, wt_stddev, wt_var, threshold_score
            ]],
            columns=[
                "sample", "read_type", "mosdepth_global_mean_cov",
                "global_mean_cov", "global_mean_cov_stddev",
                "global_mean_cov_var", f"mean_cov_at_pctile_{fmt_pctile}"
            ]
        )
        df.to_csv(output.stats, sep="\t", header=True, index=False)
    # END OF RUN BLOCK


# TODO: need wildcard comb function for
# sample and read type
rule run_all_mosdepth_assembly_read_coverage:
    input:
        check_files = expand(
            rules.mosdepth_assembly_read_coverage_window.output.check,
            sample=SAMPLES,
            read_type=["hifi", "ont"],
            aln_subset=["onlyPRI", "onlySPL"],
            mapq=MOSDEPTH_ASSM_READ_COV_MAPQ_THRESHOLDS
        )


rule run_all_mosdepth_coverage_stats:
    input:
        stats = expand(
            rules.mosdepth_coverage_stats_summary.output.stats,
            sample=SAMPLES,
            read_type=["hifi", "ont"]
        )
