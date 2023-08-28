
rule align_reads_to_complete_assembly:
    """NB: the input assembly is tagged,
    i.e. contigs carry an identifying suffix
    """
    input:
        assembly = rules.merge_and_tag_asm_units.output.mrg_fasta,
        reads = lambda wildcards: SAMPLE_INFOS[wildcards.sample][("reads", wildcards.read_type, wildcards.path_id)]
    output:
        mapped_bam = DIR_PROC.joinpath(
            "20-read-align", "00_per_file", "{sample}.{read_type}",
            "{sample}.{read_type}.{path_id}.sort.bam"
        ),
        mapped_bai = DIR_PROC.joinpath(
            "20-read-align", "00_per_file", "{sample}.{read_type}",
            "{sample}.{read_type}.{path_id}.sort.bam.bai"
        ),
        excluded_bam = DIR_PROC.joinpath(
            "20-read-align", "00_per_file", "{sample}.{read_type}",
            "{sample}.{read_type}.{path_id}.excluded.bam"
        )
    benchmark:
        DIR_RSRC.joinpath(
            "20-read-align", "00_per_file", "{sample}.{read_type}",
            "{sample}.{read_type}.{path_id}.align.rsrc"
        )
    conda:
        DIR_ENVS.joinpath("aligner", "minimap.yaml")
    threads: CPU_HIGH
    resources:
        mem_mb = lambda wildcards, attempt: 65536 + 32768 * attempt,
        time_hrs = lambda wildcards, attempt: 11 * attempt,
        sort_mem_mb = lambda wildcards, attempt: 2048 + 2048 * attempt
    params:
        readgroup = lambda wildcards: (
            f'"@RG\\tID:{wildcards.sample}_{wildcards.read_type}'
            f'\\tSM:{wildcards.sample}"'
        ),
        sam_flag_out = 1540,  # unmap, PCR-dup, QC-fail --- keep 2nd align!
        sam_threads = CPU_LOW,
        preset = lambda wildcards: f"map-{wildcards.read_type}"
    shell:
        "minimap2 -a -x {params.preset} --MD --cs --eqx -t {threads} "
        "-R {params.readgroup} {input.assembly} {input.reads}"
            " | "
        "samtools view -u -h --output-unselected {output.excluded_bam} "
        "-F {params.sam_flag_out} --threads {params.sam_threads}"
            " | "
        "samtools sort -l 9 -m {resources.sort_mem_mb}M "
        "--threads {params.sam_threads} "
        "-T {wildcards.sample}_{wildcards.read_type}_{wildcards.path_id}_mm2 "
        "-o {output.mapped_bam} "
            " && "
        "samtools index -@ {threads} {output.mapped_bam}"


# dummy rule / trigger only
localrules: collect_all_read_to_assembly_alignments
rule collect_all_read_to_assembly_alignments:
    input:
        bams = lambda wildcards: expand(
            rules.align_reads_to_complete_assembly.output.mapped_bam,
            path_id=SAMPLE_INFOS[wildcards.sample][("reads", wildcards.read_type, "path_id_all")],
            allow_missing=True
        )
    output:
        ok = DIR_PROC.joinpath("20-read-align", "00_per_file", "{sample}.{read_type}.ok")
    shell:
        "touch {output.ok}"


rule run_all_read_to_assembly_alignments:
    input:
        checkfiles = expand(
            rules.collect_all_read_to_assembly_alignments.output.ok,
            sample=SAMPLES,
            read_type=["ont", "hifi"]
        )
