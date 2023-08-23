
rule mosdepth_assembly_reference_coverage_window:
    input:
        bam = rules.minimap_assembly_to_reference_align_bam.output.bam,
        bai = rules.minimap_assembly_to_reference_align_bam.output.bai
    output:
        check = DIR_PROC.joinpath(
            "50-postprocess", "asm_ctg_refcov", "mosdepth",
            "{ref}", "{sample}.{seq_type}.mq{mapq}.wd",
            "{sample}.{seq_type}.{ref}.mq{mapq}.ok"
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


rule mosdepth_merge_region_coverage_mapq_thresholds:
    input:
        check = expand(
            rules.mosdepth_assembly_reference_coverage_window.output.check,
            mapq=MOSDEPTH_ASSM_REF_COV_MAPQ_THRESHOLDS,
            seq_type=ASSEMBLY_UNITS_NO_CONTAM,
            allow_missing=True
        )
    output:
        merged_regions = DIR_PROC.joinpath(
            "50-postprocess", "asm_ctg_refcov", "mosdepth",
            "{ref}", "{sample}.{ref}.win-ctg-cov.tsv.gz"
        )
    resources:
        mem_mb=lambda wildcards, attempt: 1024 * attempt
    run:
        import pathlib as pl
        import pandas as pd

        def assign_coverage_label(cov_value):
            if cov_value == 1:
                return "one"
            if cov_value == 0:
                return "zero"
            if cov_value < 1:
                return "low"
            if cov_value == 2:
                return "two"
            return "high"

        merged = []
        for check_file in input.check:
            regions_file = pl.Path(check_file).with_suffix(".regions.bed.gz")
            assert regions_file.is_file(), f"File not found: {regions_file}"
            name_components = regions_file.name.split(".")
            mapq_t = name_components[-4]
            assert mapq_t.startswith("mq")
            mapq_t = int(mapq_t.strip("mq"))
            asm_unit = name_components[-6]
            assert asm_unit.startswith("asm")
            regions = pd.read_csv(
                regions_file, sep="\t", header=None,
                names=["chrom", "start", "end", "ctg_align_cov"],
                index_col=["chrom", "start", "end"]
            )
            regions["label_cov"] = regions["ctg_align_cov"].apply(assign_coverage_label)
            regions.columns = pd.MultiIndex.from_tuples(
                [(asm_unit, mapq_t, c) for c in regions.columns],
                names=["asm_unit", "mapq", "statistic"]
            )
            merged.append(regions)
        if len(merged) == 1:
            merged = merged[0]
        else:
            merged = pd.concat(merged, axis=1, ignore_index=False)

        merged.to_csv(output.merged_regions, index=True, header=True, sep="\t")
    # END OF RUN BLOCK


rule run_assembly_reference_coverage:
    input:
        windowed = expand(
            rules.mosdepth_assembly_reference_coverage_window.output.check,
            ref=WILDCARDS_REF_GENOMES,
            sample=SAMPLES,
            seq_type=[f"asm-{asm_unit}" for asm_unit in ["hap1", "hap2", "unassigned", "disconnected"]],
            mapq=MOSDEPTH_ASSM_REF_COV_MAPQ_THRESHOLDS
        ),
        rdna_win = expand(
            rules.mosdepth_assembly_reference_coverage_window.output.check,
            ref=["t2tv2"],
            sample=SAMPLES,
            seq_type=["asm-rdna"],
            mapq=MOSDEPTH_ASSM_REF_COV_MAPQ_THRESHOLDS
        )

