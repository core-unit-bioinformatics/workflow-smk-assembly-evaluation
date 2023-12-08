

rule deepvariant_read_assm_alignments:
    input:
        assm = rules.merge_and_tag_asm_units.output.mrg_fasta,
        assm_idx = rules.index_merged_tagged_assembly_fasta.output.fai,
        bam = rules.merge_read_to_assembly_subset_alignments.output.bam,
        bai = rules.merge_read_to_assembly_subset_alignments.output.bai,
        clean_regions = rules.define_clean_assembly_regions.output.tag_tig
    output:
        vcfgz = DIR_PROC.joinpath(
            "60-flagging", "mismatches",
            "{sample}.{read_type}.{aln_subset}.dv-wg.vcf.gz"
        )
    benchmark:
        DIR_RSRC.joinpath(
            "60-flagging", "mismatches",
            "{sample}.{read_type}.{aln_subset}.dv-wg.rsrc"
        )
    log:
        DIR_LOG.joinpath(
            "60-flagging", "mismatches",
            "{sample}.{read_type}.{aln_subset}.dv-wg.log"
        )
    container:
        f"{CONTAINER_STORE}/{config['deepvariant']}"
    threads: CPU_MEDIUM
    resources:
        mem_mb = lambda wildcards, attempt: 24576 + 24576 * attempt,
        time_hrs = lambda wildcards, attempt: 11 * attempt,
        arch=":arch=skylake"  # docker default built with AVX512
    params:
        tempdir = lambda wildcards: DIR_PROC.joinpath(
            "temp", "deepvariant", wildcards.sample,
            wildcards.read_type, wildcards.aln_subset
        ),
        model = lambda wildcards: config["deepvariant_models"][wildcards.read_type]
    shell:
        "rm -rf {params.tempdir}"
            " && "
        "mkdir -p {params.tempdir}"
            " && "
        "/opt/deepvariant/bin/run_deepvariant --model_type {params.model} "
        "--ref {input.assm} --reads {input.bam} --num_shards {threads} "
        "--output_vcf {output.vcfgz} --regions {input.clean_regions} "
        "--noruntime_report --novcf_stats_report --sample_name {wildcards.sample} "
        "--intermediate_results_dir {params.tempdir} &> {log}"
            " ; "
        "rm -rfd {params.tempdir}"


rule apply_basic_quality_filter:
    input:
        vcf = rules.deepvariant_read_assm_alignments.output.vcfgz,
        tag_to_untag = rules.create_sed_replacement_files.output.tag_to_untag
    output:
        vcf = DIR_RES.joinpath(
            "regions", "{sample}",
            "{sample}.asmerr-{read_type}-{aln_subset}.dv-wg.vcf.gz"
        ),
        tbi = DIR_RES.joinpath(
            "regions", "{sample}",
            "{sample}.asmerr-{read_type}-{aln_subset}.dv-wg.vcf.gz.tbi"
        )
    conda:
        DIR_ENVS.joinpath("biotools", "vcftools.yaml")
    shell:
        "bcftools view -f PASS -i 'FORMAT/DP>=5' --output-type v {input.vcf} "
            " | "
        "sed -f {input.tag_to_untag}"
            " | "
        "bgzip > {output.vcf}"
            " && "
        "tabix --force --preset vcf {output.vcf}"


# rule compute_qv_estimate_mismatches:
#     input:
#         vcf = rules.apply_basic_quality_filter.output.vcf,
#         tbi = rules.apply_basic_quality_filter.output.tbi
#     output:


rule run_all_deepvariant_hifi_mismatches:
    input:
        vcf = expand(
            rules.apply_basic_quality_filter.output.vcf,
            sample=SAMPLES,
            read_type=["hifi"],
            aln_subset=["onlyPRI"]
        )
