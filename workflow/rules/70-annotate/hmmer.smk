
############### DEBUG
# Inflated HMMER/Singularity container memory footprint likely due to PBS misconfig
# This was the original scaling behavior
#        mem_mb = lambda wildcards, attempt: (65536 + 32768 * attempt) * hmmer_scaling("mem", wildcards.motif),
#        time_hrs = lambda wildcards, attempt: attempt * attempt * hmmer_scaling("time", wildcards.motif)
#######################

rule hmmer_motif_search:
    """NB: the reported hits can EITHER
    be thresholded on the E-value [-E] OR
    on the score [-T], but not both in the same run.
    The implementation here only thresholds on
    the E-value (if specified for the motif) and
    then later labels hits above the score threshold
    (if specified for the motif) as high-quality
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
    container:
        f"{CONTAINER_STORE}/{config['hmmer']}"
    threads: lambda wildcards: min(CPU_MAX, CPU_MEDIUM * hmmer_scaling("cpu", wildcards.motif))
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 102400,
        time_hrs = lambda wildcards, attempt: attempt * 48
    params:
        evalue_t = lambda wildcards: hmmer_threshold_value("evalue_t", wildcards.motif),
        nhmmer_exec = config.get("nhmmer_exec", "nhmmer")
    shell:
        "{params.nhmmer_exec} --cpu {threads} --dna "
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
