
rule hmmer_motif_search:
    """NB: the reported hits can EITHER
    be thresholded on the E-value [-E] OR
    on the score [-T], but not both in the same run.
    The implementation here only thresholds on
    the E-value (if specified for the motif) and
    then later labels hits above the score threshold
    (if specified for the motif) as high-quality

    IMPORTANT:
    Only v3.4+ of HMMER has a bug fix for an invalid
    alphabet detection:
    https://github.com/EddyRivasLab/hmmer/pull/252
    """
    input:
        asm_unit = rules.create_plain_assembly_file.output.tmp_fa,
        motif = DIR_GLOBAL_REF.joinpath("{motif}.fasta")
    output:
        txt = DIR_PROC.joinpath(
            "70-annotate", "hmmer",
            "{sample}.{asm_unit}.hmmer.wd",
            "{sample}.{asm_unit}.{motif}.hmmer-out.txt"
        ),
        table = DIR_PROC.joinpath(
            "70-annotate", "hmmer",
            "{sample}.{asm_unit}.hmmer.wd",
            "{sample}.{asm_unit}.{motif}.hmmer-table.txt"
        ),
    log:
        DIR_LOG.joinpath(
            "70-annotate", "hmmer",
            "{sample}.{asm_unit}.{motif}.hmmer.log"
        )
    benchmark:
        DIR_RSRC.joinpath(
            "70-annotate", "hmmer",
            "{sample}.{asm_unit}.{motif}.hmmer.rsrc"
        )
    conda:
        DIR_ENVS.joinpath("biotools", "hmmer.yaml")
    threads: lambda wildcards: min(CPU_MAX, CPU_LOW * hmmer_scaling("cpu", wildcards.motif))
    resources:
        mem_mb = lambda wildcards, attempt: (16384 * attempt) * hmmer_scaling("mem", wildcards.motif),
        time_hrs = lambda wildcards, attempt: attempt * attempt * hmmer_scaling("time", wildcards.motif)
    params:
        evalue_t = lambda wildcards: hmmer_threshold_value("evalue_t", wildcards.motif),
    shell:
        "nhmmer --cpu {threads} --dna "
        "-o {output.txt} --tblout {output.table} "
        "-E {params.evalue_t} "
        "{input.motif} {input.asm_unit} &> {log}"


rule run_all_hmmer_jobs:
    input:
        tables = expand(
            rules.hmmer_motif_search.output.table,
            sample=SAMPLES,
            asm_unit=ASSEMBLY_UNITS_NO_CONTAM,
            motif=HMMER_MOTIF_NAMES
        )
