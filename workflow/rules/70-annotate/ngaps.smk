
rule generate_ngaps_annotation:
    input:
        asm_units = expand(
            rules.compress_clean_assembly_sequences.output.fagz,
            asm_unit=ASSEMBLY_UNITS_NO_CONTAM,
            allow_missing=True
        )
    output:
        bed = DIR_RES.joinpath(
            "regions", "{sample}",
            "{sample}.ngaps.bed"
        )
    conda:
        DIR_ENVS.joinpath("pyseq.yaml")
    resources:
        mem_mb=lambda wildcards, attempt: 1024 * attempt
    params:
        script=find_script("localize_ngaps")
    shell:
        "{params.script} --fasta-input {input.asm_units} "
        "--output {output.bed} --name {wildcards.sample}"
