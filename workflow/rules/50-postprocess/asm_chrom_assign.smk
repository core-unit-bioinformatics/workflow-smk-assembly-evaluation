
rule estimate_chromosome_assignment:
    input:
        tsv = rules.normalize_minimap_assembly_to_reference_align_paf.output.tsv
    output:
        tsv_query = DIR_RES.joinpath(
            "reports", "ref_chrom_assign",
            "{sample}.{asm_unit}.{refgenome}.chrom-assign-by-query.tsv"
        ),
        tsv_target = DIR_RES.joinpath(
            "reports", "ref_chrom_assign",
            "{sample}.{asm_unit}.{refgenome}.chrom-assign-by-target.tsv"
        )
    conda:
        DIR_ENVS.joinpath("pyutils.yaml")
    resources:
        mem_mb = lambda wildcards, attempt: 1024 * attempt
    params:
        script=find_script("compute_chrom_assign")
    shell:
        "{params.script} --input {input.tsv} "
            "--out-query {output.tsv_query} "
            "--out-target {output.tsv_target} "


rule run_chromosome_assignments:
    input:
        tsv_query = expand(
            rules.estimate_chromosome_assignment.output.tsv_query,
            ref=["t2tv2"],
            sample=SAMPLES,
            asm_unit=[f"asm-{asm_unit}" for asm_unit in ["hap1", "hap2", "unassigned", "disconnected"]],
        ),
        tsv_target = expand(
            rules.estimate_chromosome_assignment.output.tsv_target,
            ref=["t2tv2"],
            sample=SAMPLES,
            asm_unit=[f"asm-{asm_unit}" for asm_unit in ["hap1", "hap2", "unassigned", "disconnected"]],
        ),
        tsv_qry_rdna = expand(
            rules.estimate_chromosome_assignment.output.tsv_query,
            ref=["t2tv2"],
            sample=SAMPLES,
            asm_unit=["asm-rdna"],
        ),
        tsv_trg_rdna = expand(
            rules.estimate_chromosome_assignment.output.tsv_target,
            ref=["t2tv2"],
            sample=SAMPLES,
            asm_unit=["asm-rdna"],
        )

