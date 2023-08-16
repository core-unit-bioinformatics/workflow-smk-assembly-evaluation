
rule minimap_assembly_to_reference_align_paf:
    """NB: given the wildcard composition,
    the rules in this module are not supposed
    to be used to align sequence contaminants
    to the reference (which seems like an
    illogical thing to do anyway).
    """
    input:
        ref = lambda wildcards: get_reference_assembly(wildcards.sample, wildcards.ref),
        assm = rules.compress_clean_assembly_sequences.output.fagz
    output:
        paf = DIR_RES.joinpath(
            "alignments", "contig_to_ref", "{ref}",
            "paf", "{sample}.asm-{asm_type}.{ref}.paf.gz"
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
        "pigz -p {threads} > {output.paf}"


rule normalize_minimap_assembly_to_reference_align_paf:
    input:
        paf = rules.minimap_assembly_to_reference_align_paf.output.paf
    output:
        tsv = DIR_RES.joinpath(
            "alignments", "contig_to_ref", "{ref}",
            "table", "{sample}.asm-{asm_type}.{ref}.norm-paf.tsv.gz"
        )
    conda:
        DIR_ENVS.joinpath("pyutils.yaml")
    resources:
        mem_mb = lambda wildcards, attempt: 1024 * attempt
    params:
        script=find_script("normalize_paf")
    shell:
        "{params.script} --input {input.paf} --output {output.tsv}"


rule minimap_assembly_to_reference_align_bam:
    input:
        ref = lambda wildcards: get_reference_assembly(wildcards.sample, wildcards.ref),
        assm = assm = rules.compress_clean_assembly_sequences.output.fagz
    output:
        bam = DIR_RES.joinpath(
            "alignments", "contig_to_ref", "{ref}",
            "bam", "{sample}.asm-{asm_type}.{ref}.sort.bam"
        )
        bai = DIR_RES.joinpath(
            "alignments", "contig_to_ref", "{ref}",
            "bam", "{sample}.asm-{asm_type}.{ref}.sort.bam.bai"
        )
        excluded = DIR_RES.joinpath(
            "alignments", "contig_to_ref", "{ref}",
            "bam", "{sample}.asm-{asm_type}.{ref}.unmapped.bam"
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


# TODO need general strategy to identify assembly units to use here
rule run_minimap_contig_to_ref_alignments:
    input:
        bams = expand(
                rules.minimap_assembly_to_reference_align_bam.output.bam,
                ref=WILDCARDS_REF_GENOMES,
                sample=SAMPLES,
                asm_type=["hap1", "hap2", "unassigned", "disconnected", "rdna"]
        ),
        paf = expand(
                rules.normalize_minimap_assembly_to_reference_align_paf.output.tsv,
                ref=WILDCARDS_REF_GENOMES,
                sample=SAMPLES,
                asm_type=["hap1", "hap2", "unassigned", "disconnected", "rdna"]
        )
