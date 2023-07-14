
rule minimap_assembly_to_reference_align_paf:
    input:
        ref = lambda wildcards: get_reference_assembly(wildcards.sample, wildcards.ref),
        assm = lambda wildcards: SAMPLE_INFOS[wildcards.sample][("asm", wildcards.asm_type, None)],
    output:
        paf = DIR_PROC.joinpath(
            "10-asm-align/{ref}/{sample}.asm-{asm_type}.{ref}.paf.gz"
        )
    conda:
        DIR_ENVS.joinpath("aligner", "minimap.yaml")
    threads: CPU_MEDIUM
    resources:
        time_hrs = lambda wildcards, attempt: 1 * attempt,
        mem_mb = lambda wildcards, attempt: 32768 + 32768 * attempt,
    shell:
        "minimap2 -c -x asm20 --cs --eqx -t {threads} {input.ref} {input.assm}"
            " | "
        "pigz -p > {output.paf}"


rule minimap_assembly_to_reference_align_bam:
    input:
        ref = lambda wildcards: get_reference_assembly(wildcards.sample, wildcards.ref),
        assm = lambda wildcards: SAMPLE_INFOS[wildcards.sample][("asm", wildcards.asm_type, None)],
    output:
        bam = DIR_PROC.joinpath(
            "10-asm-align/{ref}/{sample}.asm-{asm_type}.{ref}.sort.bam"
        ),
        bai = DIR_PROC.joinpath(
            "10-asm-align/{ref}/{sample}.asm-{asm_type}.{ref}.sort.bam.bai"
        ),
        exclude = DIR_PROC.joinpath(
            "10-asm-align/{ref}/{sample}.asm-{asm_type}.{ref}.excl.bam"
        )
    conda:
        DIR_ENVS.joinpath("aligner", "minimap.yaml")
    threads: CPU_MEDIUM
    resources:
        time_hrs = lambda wildcards, attempt: 1 * attempt,
        mem_mb = lambda wildcards, attempt: 32768 + 32768 * attempt,
        sort_mem_mb = lambda wildcards, attempt: 1024 * attempt
    params:
        readgroup = lambda wildcards: (
            f'"@RG\\tID:{wildcards.sample}_{wildcards.asm_type}'
            f'\\tSM:{wildcards.sample}"'
        ),
        sam_flag_out = 1540,  # unmap, PCR-dup, QC-fail --- keep 2nd align!
        sam_threads = CPU_LOW,
    shell:
        "minimap2 -a -x asm20 --cs --eqx -t {threads}"
        " -R {params.readgroup} {input.ref} {input.assm}"
            " | "
        " samtools view -u -h --output-unselected {output.exclude} "
        " -F {params.sam_flag_out} --threads {params.sam_threads}"
            " | "
        " samtools sort -l 9 -m {resources.sort_mem_mb}M "
        " --threads {params.sam_threads} "
        " -T {wildcards.sample}_{wildcards.asm_type}_{wildcards.ref}_mm2 -o {output.bam} "
            " && "
        "samtools index -@ {threads} {output.bam}"


rule run_minimap_contig_to_ref_alignments:
    input:
        bams = expand(
            DIR_PROC.joinpath(
                "10-asm-align/{ref}/{sample}.asm-{asm_type}.{ref}.sort.bam"),
                ref=WILDCARDS_REF_GENOMES,
                sample=SAMPLES,
                asm_type=["hap1", "hap2", "unassigned"]
        ),
        paf = expand(
            DIR_PROC.joinpath(
                "10-asm-align/{ref}/{sample}.asm-{asm_type}.{ref}.paf.gz"),
                ref=WILDCARDS_REF_GENOMES,
                sample=SAMPLES,
                asm_type=["hap1", "hap2", "unassigned"]
        )
