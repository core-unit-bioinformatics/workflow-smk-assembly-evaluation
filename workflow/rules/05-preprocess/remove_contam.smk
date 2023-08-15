
rule remove_assembly_contaminants:
    input:
        rep_adap = rules.normalize_merge_ncbi_fcs_adaptor_report.output.report,
        rep_contam = rules.normalize_merge_ncbi_fcs_contamination_report.output.report,
        mrg_fasta = rules.merge_and_tag_asm_units.output.mrg_fasta
    output:
        contam = DIR_PROC.joinpath(
            "05-preprocess", "remove_contam",
            "{sample}.contaminants.tmp.fa"
        ),
        asm_units = expand(
            DIR_PROC.joinpath("05-preprocess", "remove_contam", "{sample}.{asm_unit}.tmp.fa"),
            asm_unit=lambda wildcards: SAMPLE_INFOS[wildcards.sample][("asm", "all", "names")],
            allow_missing=True
        )
    log:
        DIR_LOG.joinpath("05-preprocess", "remove_contam", "{sample}.filter.log")
    resources:
        mem_mb = lambda wildcards, attempt: 2048 * attempt
    params:
        script=find_script("filter_contaminants"),
        out_pattern = lambda wildcards, output: output.contam.replace("contaminants", "asm-SEQTAG")
    shell:
        "{params.script} --input {input.mrg_fasta} --strip-tags --report "
            "--adapter-table {input.rep_adap} --contamination-table {input.rep_contam} "
            "--out-contaminants {output.contam} --out-pattern {params.out_pattern} &> {log}"


rule compress_clean_assembly_sequences:
    input:
        fasta = DIR_PROC.joinpath(
            "05-preprocess", "remove_contam",
            "{sample}.{seq_type}.tmp.fa"
        )
    output:
        fagz = DIR_RES.joinpath(
            "assemblies", "{sample}", "{sample}.{seq_type}.fasta.gz"
        ),
        fai = DIR_RES.joinpath(
            "assemblies", "{sample}", "{sample}.{seq_type}.fasta.gz.fai"
        ),
        gzi = DIR_RES.joinpath(
            "assemblies", "{sample}", "{sample}.{seq_type}.fasta.gz.gzi"
        ),
    conda:
        DIR_ENVS.joinpath("biotools.yaml")
    threads: CPU_LOW
    resources:
        mem_mb=lambda wildcards, attempt: 2048 * attempt,
    shell:
        "bgzip -c -l 6 -@ {threads} {input.fasta} > {output.fagz}"
            " && "
        "samtools faidx {output.fagz}"


# TODO fails if missing assembly unit for sample
rule run_all_clean_assembly:
    input:
        asm_units = expand(
            rules.compress_clean_assembly_sequences.output.fagz,
            sample=SAMPLES,
            seq_type=["hap1", "hap2", "unassigned", "disconnected", "ebv", "rdna", "mito"]
        ),
        contam = expand(
            rules.compress_clean_assembly_sequences.output.fagz,
            sample=SAMPLES,
            seq_type=["contaminants"]
        )
