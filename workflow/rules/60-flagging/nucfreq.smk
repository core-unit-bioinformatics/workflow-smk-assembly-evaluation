
rule build_nucfreq_cache:
    input:
        script = rules.clone_nucfreqtwo_repo.output.nucfreqtwo,
        bam = rules.merge_read_to_assembly_subset_alignments.output.bam,
        bai = rules.merge_read_to_assembly_subset_alignments.output.bai
    output:
        hdf = DIR_PROC.joinpath(
            "60-flagging", "nucfreq",
            "{sample}.{read_type}.{aln_subset}.nfdata.hdf"
        )
    benchmark:
        DIR_RSRC.joinpath(
            "60-flagging", "nucfreq",
            "{sample}.{read_type}.{aln_subset}.nucfreqtwo.rsrc"
        )
    conda:
        DIR_ENVS.joinpath("pyseq.yaml")
    threads: CPU_LOW
    resources:
        time_hrs = lambda wildcards, attempt: 8 * attempt,
        mem_mb = lambda wildcards, attempt: 24576 + 24576 * attempt
    shell:
        "{input.script} --infile {input.bam} --out-hdf-cache {output.hdf} "
        "--mode process --threads {threads} "
        "--flag-discordant-pct 10 "
        "--flag-discordant-abs 2 "
        "--flag-min-interval 500 "
        "--flag-num-hets 5 "
        "--flag-store-window 50000 "
        "--min-region-size 500000"


rule run_all_nucfreq_cache:
    input:
        hdf = expand(
            rules.build_nucfreq_cache.output.hdf,
            sample=SAMPLES,
            read_type=["hifi"],
            aln_subset=["onlyPRI"]
        )

