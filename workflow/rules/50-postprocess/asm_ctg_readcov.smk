
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
        window_size = MOSDEPTH_ASSM_REF_COV_WINDOW_SIZE,
        min_mapq = lambda wildcards: int(wildcards.mapq),
        out_prefix = lambda wildcards, output: pathlib.Path(output.check).with_suffix("")
    shell:
        "mosdepth --threads {threads} --by {params.window_size} "
            "--no-per-base --flag 0 --use-median "
            "--quantize 5:8:12:15:20:25:30: --fast-mode --mapq {params.min_mapq} "
            "{params.out_prefix} {input.bam}"
            " && "
        "touch {output.check}"


rule run_all_mosdepth_assembly_read_coverage:
    input:
        check_files = expand(
            rules.mosdepth_assembly_read_coverage_window.output.check,
            sample=SAMPLES,
            read_type=["hifi", "ont"],
            aln_subset=["onlyPRI", "onlySPL"],
            mapq=MOSDEPTH_ASSM_READ_COV_MAPQ_THRESHOLDS
        )
