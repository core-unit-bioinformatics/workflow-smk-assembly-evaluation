

# rule preprocess_gene_model:
#     input:
#         cdna = lambda wildcards: get_gene_model(wildcards.genemodel),
#         gtf = lambda wildcards: get_gene_annotation(wildcards.genemodel)
#     output:
#         full_male = DIR_LOCAL_REF.joinpath(""),
#         full_female = DIR_LOCAL_REF.joinpath("")
#     run:



rule ref_completeness_genemodel:
    input:
        fasta = lambda wildcards: get_reference_genome(wildcards.refgenome),
        cdna = lambda wildcards: get_gene_model(wildcards.genemodel)
    output:
        paf = DIR_PROC.joinpath(
            "75-completeness", "asmgene",
            "{refgenome}.{genemodel}.paf.gz"
        )
    conda:
        DIR_ENVS.joinpath("aligner", "minimap.yaml")
    threads: CPU_MEDIUM
    resources:
        time_hrs = lambda wildcards, attempt: 1 * attempt,
        mem_mb = lambda wildcards, attempt: 32768 + 32768 * attempt,
    shell:
        "minimap2 -cxsplice:hq -t {threads} {input.fasta} {input.cdna}"
            " | "
        "pigz > {output.paf}"


rule asm_completeness_genemodel:
    input:
        fasta = rules.compress_clean_assembly_sequences.output.fagz,
        cdna = lambda wildcards: get_gene_model(wildcards.genemodel)
    output:
        paf = DIR_PROC.joinpath(
            "75-completeness", "asmgene",
            "assemblies", "{sample}",
            "{sample}.{asm_unit}.{genemodel}.paf.gz"
        )
    benchmark:
        DIR_RSRC.joinpath(
            "75-completeness", "asmgene",
            "assemblies",
            "{sample}.{asm_unit}.{genemodel}.mm2.rsrc"
        )
    conda:
        DIR_ENVS.joinpath("aligner", "minimap.yaml")
    threads: CPU_MEDIUM
    resources:
        time_hrs = lambda wildcards, attempt: 1 * attempt,
        mem_mb = lambda wildcards, attempt: 32768 + 32768 * attempt,
    shell:
        "minimap2 -cxsplice:hq -t {threads} {input.fasta} {input.cdna}"
            " | "
        "pigz > {output.paf}"


rule asm_completeness_asmgene_stats:
    input:
        ref = rules.ref_completeness_genemodel.output.paf,
        assm = rules.asm_completeness_genemodel.output.paf
    output:
        stats = DIR_PROC.joinpath(
            "75-completeness", "asmgene",
            "assemblies", "{sample}",
            "{sample}.{asm_unit}.{refgenome}.{genemodel}.asmgene-stats.txt"
        )
    conda:
        DIR_ENVS.joinpath("aligner", "minimap.yaml")
    threads: 1
    resources:
        time_hrs = lambda wildcards, attempt: 1 * attempt - 1,
        mem_mb = lambda wildcards, attempt: 2048 * attempt,
    shell:
        "paftools.js asmgene -e {input.ref} {input.assm} > {output.stats}"


rule run_all_asm_completeness_asmgene:
    """
    TODO: make wildcard values for refgenome / genemodel config params
    """
    input:
        asmgene_stats = expand(
            rules.asm_completeness_asmgene_stats.output.stats,
            sample=SAMPLES,
            asm_unit=ASSEMBLY_UNITS_NO_CONTAM,
            refgenome=["t2tv2"],
            genemodel=["gencodeV43"]
        )
