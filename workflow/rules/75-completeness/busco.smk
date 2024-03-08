
rule minibusco:
    input:
        asm = rules.compress_clean_assembly_sequences.output.fagz,
        busco_db = DIR_GLOBAL_REF.joinpath("busco_db", "{odb_name}", "{odb_name}.done")
    output:
        check = DIR_PROC.joinpath(
            "75-completeness", "busco", "{sample}.{asm_unit}.{odb_name}.wd"
            "{sample}.{asm_unit}.{odb_name}.ok"
        )
    log:
        DIR_LOG.joinpath("75-completeness", "busco", "{sample}.{asm_unit}.{odb_name}.compleasm.log")
    benchmark:
        DIR_RSRC.joinpath("75-completeness", "busco", "{sample}.{asm_unit}.{odb_name}.compleasm.log")
    conda:
        DIR_ENVS.joinpath("biotools", "compleasm")
    threads: CPU_LOW
    resources:
        mem_mb = lambda wildcards, attempt: 32768 * attempt,
        time_hrs = lambda wildcards, attempt: attempt * attempt
    params:
        wd=lambda wildcards, output: pathlib.Path(output.check).parent,
        dbdir=lambda wildcards, input: pathlib.Path(input.busco_db).parent
    shell:
        "compleasm run --mode busco -L {params.dbdir} -l {wildcards.odb_name} "
        "--threads {threads} -o {params.wd} -a {input.asm} &> {log}"
            " && "
        "touch {output.check}"


# TODO: make odb db name parameter
rule run_all_minibusco:
    input:
        checks = expand(
            rules.minibusco.output.check,
            sample=SAMPLES,
            asm_unit=ASSEMBLY_UNITS_MAIN,
            odb_name=["eukaryota_odb10"]
        )




