
rule merge_read_to_assembly_subset_alignments:
    input:
        bams = lambda wildcards: expand(
            rules.filter_read_alignments_to_subset.output.bam,
            path_id=SAMPLE_INFOS[wildcards.sample][("reads", wildcards.read_type, "path_id_all")],
            allow_missing=True
        ),
        bai = lambda wildcards: expand(
            rules.filter_read_alignments_to_subset.output.bai,
            path_id=SAMPLE_INFOS[wildcards.sample][("reads", wildcards.read_type, "path_id_all")],
            allow_missing=True
        )
    output:
        bam = DIR_PROC.joinpath(
            "20-read-align", "20_merge_subsets", "{sample}.{read_type}",
            "{sample}.{read_type}.{aln_subset}.sort.bam"
        ),
        bai = DIR_PROC.joinpath(
            "20-read-align", "20_merge_subsets", "{sample}.{read_type}",
            "{sample}.{read_type}.{aln_subset}.sort.bam.bai"
        ),
    wildcard_constraints:
        aln_subset = "(" + "|".join(["onlyPRI", "onlySPL", "onlySEC"]) + ")"
    conda:
        DIR_ENVS.joinpath("aligner", "minimap.yaml")
    threads: CPU_LOW
    resources:
        mem_mb = lambda wildcards, attempt: 8192 * attempt,
        time_hrs = lambda wildcards, attempt: attempt ** 3,
    shell:
        "samtools merge -r -f --threads {threads} -o {output.bam} {input.bams}"
            " && "
        "samtools index --threads {threads} {output.bam}"
