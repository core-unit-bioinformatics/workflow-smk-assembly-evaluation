
rule normalize_contig_ref_paf_alignments:
    input:
        paf = DIR_PROC.joinpath(
            "10-asm-align/{ref}/{sample}.asm-{asm_type}.{ref}.paf.gz"
        )
    output:
        tsv = DIR_RES.joinpath(
            "alignments/contig_ref/{sample}.asm-{asm_type}.{ref}.tsv.gz"
        )
    conda:
        DIR_ENVS.joinpath("pyutils.yaml")
    resources:
        mem_mb = lambda wildcards, attempt: 1024 * attempt
    params:
        script=find_script("normalize_paf")
    shell:
        "{params.script} --input {input.paf} --output {output.tsv}"


rule estimate_chromosome_assignment:
    input:
        tsv = DIR_RES.joinpath(
            "alignments/contig_ref/{sample}.asm-{asm_type}.{ref}.tsv.gz"
        )
    output:
        tsv_query = DIR_RES.joinpath(
            "alignments/contig_ref/{sample}.asm-{asm_type}.{ref}.ctg-chrom-query.tsv"
        ),
        tsv_target = DIR_RES.joinpath(
            "alignments/contig_ref/{sample}.asm-{asm_type}.{ref}.ctg-chrom-target.tsv"
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
                "alignments/contig_ref/{sample}.asm-{asm_type}.{ref}.ctg-chrom-{view}.tsv"),
            ref=WILDCARDS_REF_GENOMES,
            sample=SAMPLES,
            asm_type=["hap1", "hap2", "unassigned"],
            view=["query", "target"]
        )

