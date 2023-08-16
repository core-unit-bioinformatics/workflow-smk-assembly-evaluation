
rule estimate_chromosome_assignment:
    input:
        tsv = rules.normalize_minimap_assembly_to_reference_align_paf.output.tsv
    output:
        tsv_query = DIR_RES.joinpath(
            "reports" "ref_chrom_assign",
            "{sample}.asm-{asm_type}.{ref}.chrom-assign-by-query.tsv"
        ),
        tsv_target = DIR_RES.joinpath(
            "reports", "ref_chrom_assign",
            "{sample}.asm-{asm_type}.{ref}.chrom-assign-by-target.tsv"
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
        tsv = expand(
            DIR_RES.joinpath(
                "reports", "ref_chrom_assign",
                "{sample}.asm-{asm_type}.{ref}.chrom-assign-by-{view}.tsv"),
            ref=WILDCARDS_REF_GENOMES,
            sample=SAMPLES,
            asm_type=["hap1", "hap2", "unassigned", "disconnected"],
            view=["query", "target"]
        ),
        tsv_rdna = expand(
            DIR_RES.joinpath(
                "reports", "ref_chrom_assign",
                "{sample}.asm-{asm_type}.{ref}.chrom-assign-by-{view}.tsv"),
            ref=["t2tv2"],
            sample=SAMPLES,
            asm_type=["rdna"],
            view=["query", "target"]
        )
