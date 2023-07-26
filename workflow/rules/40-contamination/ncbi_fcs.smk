
rule ncbi_fcs_adaptor_screening:
    input:
        sif = NCBI_FCS_ADAPTOR_SIF,
        sh_script = NCBI_FCS_ADAPTOR_SCRIPT,
        fasta = lambda wildcards: SAMPLE_INFOS[wildcards.sample][("asm", wildcards.asm_type, None)],
    output:
        check = DIR_PROC.joinpath(
            "40-contamination", "ncbi_fcs", "adaptor",
            "{sample}.asm-{asm_type}.ok")
    benchmark:
        DIR_RSRC.joinpath(
            "40-contamination", "ncbi_fcs", "adaptor",
            "{sample}.asm-{asm_type}.rsrc")
    log:
        DIR_LOG.joinpath(
            "40-contamination", "ncbi_fcs", "adaptor",
            "{sample}.asm-{asm_type}.log")
    resources:
        mem_mb = lambda wildcards, attempt: 4096 * attempt
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
        db = NCBI_FCS_GX_DB,
        py_script = NCBI_FCS_GX_SCRIPT,
        fasta = lambda wildcards: SAMPLE_INFOS[wildcards.sample][("asm", wildcards.asm_type, None)],
    output:
        check = DIR_PROC.joinpath(
            "40-contamination", "ncbi_fcs", "gx_contam",
            "{sample}.asm-{asm_type}.ok")
    benchmark:
        DIR_RSRC.joinpath(
            "40-contamination", "ncbi_fcs", "gx_contam",
            "{sample}.asm-{asm_type}.rsrc")
    log:
        DIR_LOG.joinpath(
            "40-contamination", "ncbi_fcs", "gx_contam",
            "{sample}.asm-{asm_type}.log")
    conda: DIR_ENVS.joinpath("biotools", "ncbi_fcs.yaml")
    resources:
        mem_mb = lambda wildcards, attempt: int((384 + 192 * attempt) * 1024),
        time_hrs = lambda wildcards, attempt: 47 * attempt
    params:
        out_dir = lambda wildcards, output: pathlib.Path(output.check).with_suffix(".wd"),
        tax_id = NCBI_FCS_GX_TAX_ID
    shell:
        "python3 fcs.py screen genome --fasta {input.fasta} "
        "--out-dir {params.out_dir} --gx-db {input.db} "
        "--tax-id {params.tax_id} &> {log}"
            " && "
        "touch {output.check}"


# TODO
rule run_assembly_ncbi_fcs_adaptor_screening:
    input:
        check = expand(
            rules.ncbi_fcs_adaptor_screening.output.check,
            sample=SAMPLES,
            asm_type=["disconnected"]
        )


# TODO
rule run_assembly_ncbi_fcs_gx_contamination_screening:
    input:
        check = expand(
            rules.ncbi_fcs_gx_contamination_screening.output.check,
            sample=SAMPLES,
            asm_type=["disconnected"]
        )
