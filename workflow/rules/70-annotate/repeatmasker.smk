
rule repeatmasker_assembly_run:
    """
    TODO: hard-coded default for human - make parameter
    TODO: unify wildcard naming - seq_type == asm_unit
    TODO: uses default RepeatMasker library - has to be downloaded manually
    and copied into the respective Conda env --- no offline deployment possible!
    see: github.com/rmhubley/RepeatMasker/issues/233

    NB: RepeatMasker cannot process compressed files, but the
    main process does not abort, just the I/O part seems to die
    """
    input:
        fasta = rules.compress_clean_assembly_sequences.output.fagz,
    output:
        repmask_out = multiext(
            str(DIR_PROC.joinpath(
                "70-annotate", "repeatmasker",
                "{sample}.{seq_type}.repmask.wd",
                "{sample}.{seq_type}.repmask.tmp.fa"
            )),
            ".cat.gz", ".masked", ".out", ".tbl"
        )
    log:
        DIR_LOG.joinpath(
            "70-annotate", "repeatmasker",
            "{sample}.{seq_type}.repmask.log"
        )
    benchmark:
        DIR_RSRC.joinpath(
            "70-annotate", "repeatmasker",
            "{sample}.{seq_type}.repmask.rsrc"
        )
    conda:
        DIR_ENVS.joinpath("biotools", "motifs.yaml")
    threads: CPU_MEDIUM
    resources:
        mem_mb = lambda wildcards, attempt, input: attempt * get_repeatmasker_run_memory_mb(input.size_mb, compressed=True),
        time_hrs = lambda wildcards, attempt, input: attempt * get_repeatmasker_run_time_hrs(input.size_mb, compressed=True),
    params:
        out_dir = lambda wildcards, output: pathlib.Path(output.repmask_out[0]).parent,
        unzip_tmp = lambda wildcards, output: pathlib.Path(output.repmask_out[0]).parent.with_suffix(".tmp.fa"),
    shell:
        "pigz -p {threads} -d -c {input.fasta} > {params.unzip_tmp}"
            " && "
        "RepeatMasker -pa {threads} -s -dir {params.out_dir} "
        "-species human {params.unzip_tmp} &> {log}"
            " ; "
        "rm -f {params.unzip_tmp}"


rule collect_repeatmasker_output:
    input:
        repmask_files = rules.repeatmasker_assembly_run.output.repmask_out
    output:
        tar = DIR_RES.joinpath(
                "annotations", "repeatmasker",
                "{sample}.{seq_type}.repmask.raw.tar.gz",
            )
    resources:
        time_hrs=lambda wildcards, attempt: max(0, attempt - 1)
    run:
        import pathlib as pl
        import subprocess as sp

        tar_cmd = f"tar czf {output.tar} -C "
        for idx, filepath in enumerate(input.repmask_files, start=0):
            if idx == 0:
                tar_dir = pl.Path(filepath).parent
                tar_cmd += f"{tar_dir} "
            filename = pl.Path(filepath).name
            tar_cmd += f"./{filename} "
        try:
            _ = sp.check_output(tar_cmd, shell=True)
        except sp.CalledProcessError as spe:
            logerr(
                (
                    "tar of RepeatMasker output failed: "
                    f"{wildcards.sample} / {wildcards.seq_type} "
                    "\nOriginal command string:\n"
                    f"{tar_cmd}"
                ))
            raise spe
    # END OF RUN BLOCK


rule normalize_repeatmasker_table:
    input:
        table = DIR_PROC.joinpath(
            "70-annotate", "repeatmasker",
            "{sample}.{seq_type}.repmask.wd",
            "{sample}.{seq_type}.repmask.tmp.fa.out"
        ),
    output:
        tsv = DIR_RES.joinpath(
            "annotations", "repeatmasker",
            "{sample}.{seq_type}.repmask.matches.tsv.gz"
        ),
    resources:
        mem_mb = lambda wildcards, attempt: 1024 * attempt
    run:
        import gzip
        import pandas as pd

        names_explained = [
            (0, 'sw_score', 'Smith-Waterman score of the match'),
            (1, 'pct_subs', 'Pct. substitutions in matching region compared to the consensus'),
            (2, 'pct_del', 'Pct. of bases opposite a gap in the query sequence (deleted bp)'),
            (3, 'pct_ins', 'Pct. of bases opposite a gap in the repeat consensus (inserted bp)'),
            (4, 'query_name', 'Query sequence name'),
            (5, 'query_start', 'Starting position of match in query sequence'),
            (6, 'query_end', 'Ending position of match in query sequence'),
            (7, 'query_after_bp', 'bp in query sequence past the ending position of match'),
            (8, 'is_complement_match', 'True: match is with the Complement of the consensus sequence in the database'),
            (9, 'repeat_name', 'Name of the matching interspersed repeat'),
            (10, 'repeat_class', 'Class of the repeat'),
            (11, 'repeat_start', 'Starting position of match in database sequence (using top-strand numbering)'),
            (12, 'repeat_end', 'Ending position of match in database sequence'),
            (13, 'repeat_before_bp', 'bp in (complement of) the repeat consensus sequence prior to beginning of the match'),
            (14, 'match_ID', 'Consecutive ID of match'),
            (15, 'is_partly_included', 'True: there is a higher-scoring match whose domain partly (<80%) includes the domain of this match')
        ]

        df = pd.read_csv(
            input.table,
            header=None,
            index_col=False,
            delimiter=r"\s+",
            skip_blank_lines=True,
            skiprows=[0,1],
            comment='#',
            names=[n[1] for n in names_explained]
        )

        for c in ['query_after_bp', 'repeat_start', 'repeat_before_bp']:
            df[c] = df[c].apply(lambda x: int(x.strip('()')))
            df[c] = df[c].astype(int)

        df['is_complement_match'] = df['is_complement_match'].apply(lambda x: True if x == 'C' else False)
        df['is_complement_match'] = df['is_complement_match'].astype(bool)
        df['is_partly_included'] = df['is_partly_included'].apply(lambda x: True if x == '*' else False)
        df['is_partly_included'] = df['is_partly_included'].astype(bool)
        df.sort_values(['query_name', 'query_start', 'sw_score'], ascending=[True, True, False], inplace=True)

        with gzip.open(output.tsv, 'wt') as table:
            for pos, head, comment in names_explained:
                _ = table.write(f'##.{pos} - {head} - {comment}\n')
        df.to_csv(output.tsv, mode="a", sep='\t', header=True, index=False)
    # END OF RUN BLOCK


rule run_all_repeatmasker_jobs:
    """NB: failed RepeatMasker runs do not clean
    up after themselves (presumably to keep debugging
    information intact) and there is no switch to change
    that behavior. Hence, this trigger rule also performs
    the cleanup operation. Obviously, this can only be
    kicked off after all jobs have passed (restarted until
    completion), which makes this a typical candidate for
    leaving behind garbage in case the pipeline is
    interrupted in some way. Extremely annoying!!!
    """
    input:
        norm_tables = expand(
            rules.normalize_repeatmasker_table.output.tsv,
            sample=SAMPLES,
            seq_type=ASSEMBLY_UNITS_NO_CONTAM
        ),
        tars = expand(
            rules.collect_repeatmasker_output.output.tar,
            sample=SAMPLES,
            seq_type=ASSEMBLY_UNITS_NO_CONTAM
        ),
    shell:
        "rm -rf RM_*"

