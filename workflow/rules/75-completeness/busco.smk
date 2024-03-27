
localrules: compleasm_list_available_busco_lineages
rule compleasm_list_available_busco_lineages:
    """NB: this makes use of the wildcard to allow for later
    db updates w/o touching all of the already generated outputs
    """
    input:
        busco_db = DIR_GLOBAL_REF.joinpath("busco_db")
    output:
        lineage = DIR_PROC.joinpath(
            "75-completeness", "busco", "lineage_{odb_name}.check.txt"
        )
    conda:
        DIR_ENVS.joinpath("biotools", "compleasm.yaml")
    shell:
        "compleasm list --local --library_path {input.busco_db} > {output.lineage}"


localrules: confirm_busco_lineage_exists
rule confirm_busco_lineage_exists:
    input:
        lineage = rules.compleasm_list_available_busco_lineages.output.lineage
    output:
        exists = DIR_PROC.joinpath(
            "75-completeness", "busco", "lineage_{odb_name}.exists.txt"
        )
    run:
        lineage_found = False
        with open(input.lineage, "r") as listing:
            first_line = listing.readline()
            if not first_line.startswith("Local available lineages"):
                err_msg = (
                    "Unexpected first line in compleasm / busco lineage file:\n"
                    f"{input.lineage}\n"
                    f"First line: {first_line.strip()}"
                )
                logerr(err_msg)
                raise ValueError(err_msg)
            for line in listing:
                if wildcards.odb_name in line or line.strip() in wildcards.odb_name:
                    lineage_found = True
                    break
        if not lineage_found:
            err_msg = "Requested BUSCO lineage not locally available: {wildcards.odb_name}"
            logerr(err_msg)
            raise ValueError(err_msg)
        with open(output.exists, "w") as dump:
            _ = dump.write(get_timestamp() + "\n")
            _ = dump.write(wildcards.odb_name + "\n")
    # END OF RUN BLOCK


rule compleasm_busco_mode:
    input:
        asm = get_asm_unit,
        lineage_exists = rules.confirm_busco_lineage_exists.output.exists,
        busco_db = DIR_GLOBAL_REF.joinpath("busco_db")
    output:
        summary = DIR_PROC.joinpath(
            "75-completeness", "busco", "{sample}.{asm_unit}.{odb_name}.wd",
            "summary.txt"
        )
    log:
        DIR_LOG.joinpath("75-completeness", "busco", "{sample}.{asm_unit}.{odb_name}.compleasm.log")
    benchmark:
        DIR_RSRC.joinpath("75-completeness", "busco", "{sample}.{asm_unit}.{odb_name}.compleasm.log")
    conda:
        DIR_ENVS.joinpath("biotools", "compleasm.yaml")
    threads: CPU_LOW
    resources:
        mem_mb = lambda wildcards, attempt: 32768 * attempt,
        time_hrs = lambda wildcards, attempt: attempt * attempt
    params:
        wd=lambda wildcards, output: pathlib.Path(output.summary).parent,
    shell:
        "compleasm run --mode busco -L {input.busco_db} -l {wildcards.odb_name} "
        "--threads {threads} -o {params.wd} -a {input.asm} &> {log}"


# TODO: make odb db name parameter
rule run_all_compleasm:
    input:
        checks = expand(
            rules.compleasm_busco_mode.output.summary,
            sample=SAMPLES,
            asm_unit=ASSEMBLY_UNITS_MAIN,
            odb_name=["eukaryota_odb10", "primates_odb10"]
        )




