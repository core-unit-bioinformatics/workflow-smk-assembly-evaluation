
rule ncbi_fcs_adaptor_screening:
    input:
        sif = NCBI_FCS_ADAPTOR_SIF,
        sh_script = NCBI_FCS_ADAPTOR_SCRIPT,
        fasta = rules.merge_and_tag_asm_units.output.mrg_fasta
    output:
        check = DIR_PROC.joinpath(
            "05-preprocess", "ncbi_fcs",
            "{sample}.adaptor.ok")
    benchmark:
        DIR_RSRC.joinpath(
            "05-preprocess", "ncbi_fcs",
            "{sample}.adaptor.rsrc")
    log:
        DIR_LOG.joinpath(
            "05-preprocess", "ncbi_fcs",
            "{sample}.adaptor.log")
    resources:
        mem_mb = lambda wildcards, attempt: 8192 * attempt,
        time_hrs = lambda wildcards, attempt: attempt
    params:
        out_dir = lambda wildcards, output: pathlib.Path(output.check).with_suffix(".wd"),
        taxonomy = NCBI_FCS_ADAPTOR_TAXONOMY
    shell:
        "mkdir -p {params.out_dir}"
            " && "
        "{input.sh_script} --fasta-input {input.fasta} --output-dir {params.out_dir} "
        "--{params.taxonomy} --container-engine singularity --image {input.sif} "
            " && "
        "touch {output.check}"


rule ncbi_fcs_gx_contamination_screening:
    input:
        sif = NCBI_FCS_GX_SIF,
        db = NCBI_FCS_GX_DB_PATH,
        py_script = NCBI_FCS_GX_SCRIPT,
        fasta = rules.merge_and_tag_asm_units.output.mrg_fasta
    output:
        check = DIR_PROC.joinpath(
            "05-preprocess", "ncbi_fcs",
            "{sample}.gx-contam.ok")
    benchmark:
        DIR_RSRC.joinpath(
            "05-preprocess", "ncbi_fcs",
            "{sample}.gx-contam.rsrc")
    log:
        DIR_LOG.joinpath(
            "05-preprocess", "ncbi_fcs",
            "{sample}.gx_contam.log")
    conda: DIR_ENVS.joinpath("biotools", "ncbi_fcs.yaml")
    resources:
        mem_mb = lambda wildcards, attempt: int((384 + 192 * attempt) * 1024),
        time_hrs = lambda wildcards, attempt: 23 * attempt
    params:
        out_dir = lambda wildcards, output: pathlib.Path(output.check).with_suffix(".wd"),
        tax_id = NCBI_FCS_GX_TAX_ID,
        db_name = NCBI_FCS_GX_DB_NAME
    shell:
        "export FCS_DEFAULT_IMAGE={input.sif} ; "
        "python3 {input.py_script} screen genome --fasta {input.fasta} "
        "--out-dir {params.out_dir} --gx-db {input.db}/{params.db_name} "
        "--tax-id {params.tax_id} &> {log}"
            " && "
        "touch {output.check}"


# TODO
rule run_assembly_ncbi_fcs_adaptor_screening:
    input:
        check = expand(
            rules.ncbi_fcs_adaptor_screening.output.check,
            sample=SAMPLES,
        )


# TODO
rule run_assembly_ncbi_fcs_gx_contamination_screening:
    input:
        check = expand(
            rules.ncbi_fcs_gx_contamination_screening.output.check,
            sample=SAMPLES,
        )
