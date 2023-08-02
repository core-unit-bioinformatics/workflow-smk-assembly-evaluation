
rule create_sequence_tags_file:
    input:
        asm_seq = lambda wildcards: SAMPLE_INFOS[wildcards.sample][("asm", "all", "files")]
    output:
        tagfile = DIR_PROC.joinpath(
            "05-preprocess", "merge_tag_asm", "{sample}.asm-tags.tsv"
        )
    run:
        import io
        import pathlib as pl

        # by construction, these lists are ordered and matched
        tag_names = SAMPLE_INFOS[wildcards.sample][("asm", "all", "names")]
        seq_files = SAMPLE_INFOS[wildcards.sample][("asm", "all", "files")]
        assert len(tag_names) == len(seq_files)

        tags = io.StringIO()
        for seqfile, tag in zip(seq_files, tag_names):
            filename = pl.Path(seqfile).name
            tags.write(f"{filename}\t{tag}\n")

        with open(output.tagfile, "w") as dump:
            _ = dump.write(tags.getvalue())
    # END OF RUN BLOCK


rule merge_and_tag_asm_units:
    input:
        asm_seq = lambda wildcards: SAMPLE_INFOS[wildcards.sample][("asm", "all", "files")],
        tags = DIR_PROC.joinpath(
            "05-preprocess", "merge_tag_asm", "{sample}.asm-tags.tsv"
        )
    output:
        mrg_fasta = DIR_PROC.joinpath(
            "05-preprocess", "merge_tag_asm", "{sample}.asm-mrg-tag.fasta"
        )
    log:
        DIR_LOG.joinpath(
            "05-preprocess", "merge_tag_asm", "{sample}.asm-mrg-tag.log"
        )
    conda:
        DIR_ENVS.joinpath("pyutils.yaml")
    resources:
        mem_mb = lambda wildcards, attempt: 2048 * attempt,
        time_hrs = lambda wildcards, attempt: attempt * attempt
    params:
        script = find_script("fasta_tag_merge")
    shell:
        "{params.script} --input {input.asm_seq} --seq-tags {input.tags} "
            "--report --output {output.mrg_fasta} 2> "
