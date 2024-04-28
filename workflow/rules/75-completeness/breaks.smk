
rule find_contig_alignment_breaks:
    """The output in query coordinates can be
    used as is in terms of an annotation file
    (dump to results folder).
    The output in reference coordinates must (should)
    be postprocessed to merge the results for all
    samples to then create a global annotation for
    remaining gaps in the assemblies.
    """
    input:
        norm_paf = get_contig_to_reference_norm_paf
    output:
        trg_all = DIR_PROC.joinpath(
            "75-completeness", "breaks", "{refgenome}",
            "{sample}.{asm_unit}.{refgenome}.{aln_type}.ctgaln-label.trg-all.bed.gz"
        ),
        trg_cmp = DIR_PROC.joinpath(
            "75-completeness", "breaks", "{refgenome}",
            "{sample}.{asm_unit}.{refgenome}.{aln_type}.ctgaln-label.trg-breaks.bed.gz"
        ),
        qry_all = DIR_RES.joinpath(
            "regions", "{sample}",
            "{sample}.{asm_unit}.{refgenome}.{aln_type}.ctgaln-label.qry-all.bed.gz"
        ),
        qry_cmp = DIR_RES.joinpath(
            "regions", "{sample}",
            "{sample}.{asm_unit}.{refgenome}.{aln_type}.ctgaln-label.qry-breaks.bed.gz"
        ),
    wildcard_constraints:
        aln_type = "(blevel|coarse)"
    conda:
        DIR_ENVS.joinpath("pyutils.yaml")
    resources:
        mem_mb=lambda wildcards, attempt: 1024 * attempt
    params:
        script=find_script("find_discontigs"),
        query_label=lambda wildcards: f"{wildcards.sample}.{wildcards.asm_unit}"
    shell:
        "{params.script} --alignments {input.norm_paf} --dump-bed-like "
        "--target-label {wildcards.refgenome} --query-label {params.query_label} "
        "--out-target-all {output.trg_all} --out-target-complement {output.trg_cmp} "
            " && "
        "{params.script} --alignments {input.norm_paf} --dump-bed-like "
        "--target-label {wildcards.refgenome} --query-label {params.query_label} "
        "--add-label-description "
        "--out-query-all {output.qry_all} --out-query-complement {output.qry_cmp} "


rule run_all_label_contig_alignments:
    input:
        ctgaln_labeled = expand(
            rules.find_contig_alignment_breaks.output,
            sample=SAMPLES,
            asm_unit=ASSEMBLY_UNITS_MAIN,
            refgenome=WILDCARDS_REF_GENOMES,
            aln_type=["blevel", "coarse"]
        )
