

localrules: define_clean_assembly_regions
rule define_clean_assembly_regions:
    input:
        asm_tags = rules.create_sequence_tags_file.output.tagfile,
        fai_files = expand(
            rules.compress_clean_assembly_sequences.output.fai,
            seq_type=ASSEMBLY_UNITS_NO_CONTAM,
            allow_missing=True,
        )
    output:
        bed = DIR_RES.joinpath(
            "regions", "{sample}", "{sample}.uncontaminated.bed"
        )
    run:
        import pandas as pd
        import pathlib as pl

        tags = []
        with open(input.asm_tags, "r") as asm_tags:
            for line in asm_tags:
                if not line.strip():
                    continue
                filename, tag = line.strip().split()
                tags.append(tag)
        assert len(tags) == len(set(tags))

        regions = []
        for fai_file in input.fai_files:
            fai_name = pl.Path(fai_file).name
            fai_tag = [tag for tag in tags if tag in fai_name]
            if len(fai_tag) != 1:
                raise ValueError(f"Cannot tag FASTA index file: {fai_name}")
            fai_tag = fai_tag[0]
            fai_regions = pd.read_csv(
                fai_file, sep="\t", header=None,
                names=["contig", "end"], usecols=[0, 1]
            )
            fai_regions["start"] = 0
            fai_regions["contig"] += f".{fai_tag}"
            regions.append(fai_regions)

        regions = pd.concat(regions, axis=0, ignore_index=False)
        regions.sort_values(["contig", "start"], inplace=True)

        with open(output.bed) as dump:
            _ = dump.write("#")
            regions[["contig", "start", "end"]].to_csv(
                dump, sep="\t", header=True, index=False
            )
    # END OF RUN BLOCK


rule filter_read_alignments_to_subset:
    """ This is designed to keep only
    the primary alignments (for tools
    such as NucFreq) or the supplementary
    alignments for downstream mosdepth
    windowed coverage computation.
    """
    input:
        bam = rules.align_reads_to_complete_assembly.output.mapped_bam,
        bai = rules.align_reads_to_complete_assembly.output.mapped_bai,
        clean_regions = rules.define_clean_assembly_regions.output.bed
    output:
        bam = DIR_PROC.joinpath(
            "20-read-align", "10_filter_per_file", "{sample}.{read_type}",
            "{sample}.{read_type}.{path_id}.{aln_subset}.sort.bam"
        ),
        bai = DIR_PROC.joinpath(
            "20-read-align", "10_filter_per_file", "{sample}.{read_type}",
            "{sample}.{read_type}.{path_id}.{aln_subset}.sort.bam.bai"
        ),
    benchmark:
        DIR_RSRC.joinpath(
            "20-read-align", "10_filter_per_file", "{sample}.{read_type}",
            "{sample}.{read_type}.{path_id}.{aln_subset}.samtools.rsrc"
        )
    conda:
        DIR_ENVS.joinpath("aligner", "minimap.yaml")
    threads: CPU_LOW
    resources:
        mem_mb = lambda wildcards, attempt: {
            "onlySPL": 12288 + 12288 * attempt,
            "onlySEC": 12288 + 12288 * attempt,
            "onlyPRI": 24576 + 24576 * attempt
        }[wildcards.aln_subset],
        time_hrs = lambda wildcards, attempt: attempt * attempt,
        sort_mem_mb = lambda wildcards, attempt: 2048 + 2048 * attempt
    params:
        select_flag = lambda wildcards: {
            "onlyPRI": "--exclude-flags 3844",
            "onlySPL": "--require-flags 2048",
            "onlySEC": "--require-flags 256",
        }[wildcards.aln_subset],
        temp_prefix = lambda wildcards: f"temp/20-read-align/{wildcards.path_id}"
        # onlyPRI: unmapped; not primary; QC fail; duplicate; supplementary
        # onlySPL: supplementary
        # onlySEC: not primary
    shell:
        "mkdir -p {params.temp_prefix}"
            " && "
        "samtools view --threads {threads} -u {params.select_flag} "
        "--region-file {input.clean_regions} {input.bam} "
            " | "
        "samtools sort -l 9 -m {resources.sort_mem_mb}M --threads {threads} "
        "-T {params.temp_prefix}/{wildcards.sample}_{wildcards.read_type}_{wildcards.aln_subset}_mm2 "
        "-o {output.bam}"
            " && "
        "samtools index -@ {threads} {output.bam}"
            " ; "
        "rm -rfd {params.temp_prefix}"

