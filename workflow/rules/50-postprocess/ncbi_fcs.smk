
rule normalize_merge_ncbi_fcs_adaptor_report:
    input:
        check = DIR_PROC.joinpath(
            "40-contamination", "ncbi_fcs", "adaptor",
            "{sample}.asm-{asm_type}.ok"
        ),
        fasta = lambda wildcards: SAMPLE_INFOS[wildcards.sample][("asm", wildcards.asm_type, None)],
    output:
        report = DIR_RES.joinpath(
            "contamination", "{sample}.asm-{asm_type}.fcs-report-adaptor.tsv"
        ),
        stats = DIR_RES.joinpath(
            "contamination", "{sample}.asm-{asm_type}.fcs-report-adaptor.stats.tsv"
        )
    conda: DIR_ENVS.joinpath("pyutils.yaml")
    resources:
        mem_mb = lambda wildcards, attempt: 1024 * attempt,
        time_hrs = lambda wildcards, attempt: attempt
    params:
        rep_file = lambda wildcards, input: pathlib.Path(
            input.check).with_suffix(".wd").joinpath("fcs_adaptor_report.txt"),
        script = find_script("normalize_merge_report.py")
    shell:
        "{params.script} --adaptor --report {params.rep_file} "
            "--fasta {input.fasta} --table {output.report} "
            "--statistics {output.stats}"



def load_contam_report_file(check_file_path):
    """TODO remove this
    """

    fcs_wd = pathlib.Path(check_file_path).with_suffix(".wd")
    if not fcs_wd.is_dir():
        report_file = fcs_wd.joinpath("fcs-run-incomplete.txt")
    else:
        wd_files = list(fcs_wd.glob("*.fcs_gx_report.txt"))
        assert len(wd_files) == 1
        report_file = wd_files[0]
    return report_file


rule normalize_merge_ncbi_fcs_contamination_report:
    input:
        check = DIR_PROC.joinpath(
            "40-contamination", "ncbi_fcs", "gx_contam",
            "{sample}.asm-{asm_type}.ok"
        ),
        fasta = lambda wildcards: SAMPLE_INFOS[wildcards.sample][("asm", wildcards.asm_type, None)],
    output:
        report = DIR_RES.joinpath(
            "contamination", "{sample}.asm-{asm_type}.fcs-report-contam.norm.tsv"
        ),
        stats = DIR_RES.joinpath(
            "contamination", "{sample}.asm-{asm_type}.fcs-report-contam.stats.tsv"
        )
    conda: DIR_ENVS.joinpath("pyutils.yaml")
    resources:
        mem_mb = lambda wildcards, attempt: 1024 * attempt,
        time_hrs = lambda wildcards, attempt: attempt
    params:
        rep_file = lambda wildcards, input: load_contam_report_file(input.check)
    shell:
        "{params.script} --contamination --report {params.rep_file} "
            "--fasta {input.fasta} --table {output.report} "
            "--statistics {output.stats}"


rule run_all_ncbi_fcs_reports:
    input:
        reports = expand(
            DIR_RES.joinpath(
            "contamination", "{sample}.asm-{asm_type}.fcs-report-{report_type}.norm.tsv"),
            sample=SAMPLES,
            asm_type=["disconnected", "unassigned"],
            report_type=["contam", "adaptor"]
        )
