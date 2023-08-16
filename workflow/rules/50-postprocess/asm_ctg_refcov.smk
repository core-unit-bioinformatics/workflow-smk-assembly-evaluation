
rule mosdepth_assembly_reference_coverage_window:
    input:
        bam = rules.minimap_assembly_to_reference_align_bam.output.bam,
        bai = rules.minimap_assembly_to_reference_align_bam.output.bai
    output:
        check = DIR_PROC.joinpath(
            "50-postprocess", "asm_ctg_refcov", "mosdepth",
            "{ref}", "{sample}.asm-{asm_type}.mq{mapq}.wd",
            "{sample}.asm-{asm_type}.{ref}.mq{mapq}.ok"
        )
    threads: CPU_LOW
    conda: DIR_ENVS.joinpath("biotools", "mosdepth.yaml")
    resources:
        mem_mb = lambda wildcards, attempt: 4096 * attempt,
        time_hrs = lambda wildcards, attempt: max(0, attempt-1)
    params:
        window_size = MOSDEPTH_ASSM_REF_COV_WINDOW_SIZE,
        exclude_flag = 1796,  # same as mosdepth default
        min_mapq = lambda wildcards: int(wildcards.mapq),
        out_prefix = lambda wildcards, output: pathlib.Path(output.check).with_suffix("")
    shell:
        "mosdepth --threads {threads} --by {params.window_size} "
            "--no-per-base --flag {params.exclude_flag} --use-median "
            "--quantize 0:1:2:3:4:5:6: --fast-mode --mapq {params.min_mapq} "
            "{params.out_prefix} {input.bam}"
            " && "
        "touch {output.check}"



rule run_assembly_reference_coverage:
    input:
        windowed = expand(
            rules.mosdepth_assembly_reference_coverage_window.output.check,
            ref=WILDCARDS_REF_GENOMES,
            sample=SAMPLES,
            asm_type=["hap1", "hap2", "unassigned", "disconnected"],
            mapq=["00", "60"]
        ),
        rdna_win = expand(
            rules.mosdepth_assembly_reference_coverage_window.output.check,
            ref=["t2tv2"],
            sample=SAMPLES,
            asm_type=["rdna"],
            mapq=["00", "60"]
        )

