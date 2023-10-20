
rule repeatmasker_assembly_run:
    """
    TODO: hard-coded default for human - make parameter
    TODO: unify wildcard naming - seq_type == asm_unit
    TODO: uses default RepeatMasker library - has to be downloaded manually
    and copied into the respective Conda env --- no offline deployment possible!

    NB: RepeatMasker cannot process compressed files, but the
    main process does not abort, just the I/O part seems to die
    """
    input:
        fasta = rules.compress_clean_assembly_sequences.output.fagz,
    output:
        check = DIR_PROC.joinpath(
            "70-annotate", "repeatmasker",
            "{sample}.{seq_type}.repmask.ok"
        )
        # repmask_output = multiext(
        #     'output/motif_search/40_repmask/{sample}.{hifi_type}.{ont_type}.{mapq}.chrY/{sample}.{hifi_type}.{ont_type}.{mapq}.chrY',
        #     '.fasta.cat.gz', '.fasta.masked', '.fasta.out', '.fasta.tbl'
        # )
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
        time_hrs = lambda wildcards, attempt, input: attempt * get_repeatmasker_run_memory_mb(input.size_mb, compressed=True),
    params:
        out_dir = lambda wildcards, output: pathlib.Path(output.check).with_suffix(".wd"),
        unzip_tmp = lambda wildcards, output: pathlib.Path(output.check).with_suffix(".tmp.fa"),
    shell:
        "pigz -p {threads} -d -c {input.fasta} > {params.unzip_tmp}"
            " && "
        "RepeatMasker -pa {threads} -s -dir {params.out_dir} "
        "-species human {params.unzip_tmp} &> {log}"
            " && "
        "touch {output.check}"
            " ; "
        "rm -f {params.unzip_tmp}"


rule run_all_repeatmasker_jobs:
    input:
        check_files = expand(
            rules.repeatmasker_assembly_run.output.check,
            sample=SAMPLES,
            seq_type=ASSEMBLY_UNITS_NO_CONTAM
        )


### BELOW
# to be adapted - copied from chrY pipeline

# rule repmsk_ref_chry_contigs:
#     """
#     Same as for HMMER above.
#     For reference chrY, no renaming needed, that part is just kept here
#     to not break anything else
#     """
#     input:
#         fasta = 'references_derived/{sample}_chrY.fasta',
#         rename = 'output/subset_wg/25_name_mappings/{sample}.{hifi_type}.{ont_type}.na.chrY.names.nto-map.sed'
#     output:
#         tmp_fasta = temp('output/motif_search/40_repmask/{sample}.{hifi_type}.{ont_type}.na.chrY.fasta'),
#         repmask_output = multiext(
#             'output/motif_search/40_repmask/{sample}.{hifi_type}.{ont_type}.na.chrY/{sample}.{hifi_type}.{ont_type}.na.chrY',
#             '.fasta.cat.gz', '.fasta.masked', '.fasta.out', '.fasta.tbl'
#         )
#     log:
#         'log/output/motif_search/40_repmask/{sample}.{hifi_type}.{ont_type}.na.chrY.repmask.log',
#     benchmark:
#         'rsrc/output/motif_search/40_repmask/{sample}.{hifi_type}.{ont_type}.na.chrY.repmask.rsrc',
#     conda:
#         '../envs/biotools.yaml'
#     wildcard_constraints:
#         sample = '(T2T|GRCh38)'
#     threads: config['num_cpu_low']
#     resources:
#         mem_mb = lambda wildcards, attempt: 12288 * attempt,
#         walltime = lambda wildcards, attempt: f'{attempt*attempt:02}:59:00',
#     params:
#         out_dir = 'output/motif_search/40_repmask/{sample}.{hifi_type}.{ont_type}.na.chrY'
#     shell:
#         'sed -f {input.rename} {input.fasta} > {output.tmp_fasta} && '
#         'RepeatMasker -pa {threads} -s -dir {params.out_dir} -species human {output.tmp_fasta} &> {log}'


# rule collect_repeatmasker_output:
#     input:
#         report = 'output/motif_search/40_repmask/{sample}.{hifi_type}.{ont_type}.{mapq}.chrY/{sample}.{hifi_type}.{ont_type}.{mapq}.chrY.fasta.tbl',
#         masked = 'output/motif_search/40_repmask/{sample}.{hifi_type}.{ont_type}.{mapq}.chrY/{sample}.{hifi_type}.{ont_type}.{mapq}.chrY.fasta.masked',
#         aln = 'output/motif_search/40_repmask/{sample}.{hifi_type}.{ont_type}.{mapq}.chrY/{sample}.{hifi_type}.{ont_type}.{mapq}.chrY.fasta.cat.gz',
#         table = 'output/motif_search/40_repmask/{sample}.{hifi_type}.{ont_type}.{mapq}.chrY/{sample}.{hifi_type}.{ont_type}.{mapq}.chrY.fasta.out',
#         rename = 'output/subset_wg/25_name_mappings/{sample}.{hifi_type}.{ont_type}.{mapq}.chrY.names.otn-map.sed',
#     output:
#         fasta_rn = 'output/motif_search/45_rm_norm/{sample}/{sample}.{hifi_type}.{ont_type}.{mapq}.chrY.rm-mask.fasta',
#         result_tar = 'output/motif_search/45_rm_norm/{sample}/{sample}.{hifi_type}.{ont_type}.{mapq}.chrY.rm-out.tar.gz',
#     conda:
#         '../envs/biotools.yaml'
#     resources:
#         mem_mb = lambda wildcards, attempt: 1024 * attempt
#     params:
#         tar_dir = lambda wildcards, input: pathlib.Path(input.report).parent,
#         tar_report = lambda wildcards, input: pathlib.Path(input.report).name,
#         tar_masked = lambda wildcards, input: pathlib.Path(input.masked).name,
#         tar_aln = lambda wildcards, input: pathlib.Path(input.aln).name,
#         tar_table = lambda wildcards, input: pathlib.Path(input.table).name,
#     shell:
#         'seqtk seq -C -A {input.masked} | sed -f {input.rename} > {output.fasta_rn}'
#             ' && '
#         'tar czf {output.result_tar} -C {params.tar_dir} ./{params.tar_report} ./{params.tar_masked} '
#             './{params.tar_aln} ./{params.tar_table}'


# rule normalize_repeatmasker_table:
#     input:
#         table = 'output/motif_search/40_repmask/{sample}.{hifi_type}.{ont_type}.{mapq}.chrY/{sample}.{hifi_type}.{ont_type}.{mapq}.chrY.fasta.out',
#         rename = 'output/subset_wg/25_name_mappings/{sample}.{hifi_type}.{ont_type}.{mapq}.chrY.names.otn-map.json',
#     output:
#         tsv = 'output/motif_search/45_rm_norm/{sample}/{sample}.{hifi_type}.{ont_type}.{mapq}.chrY.matches.tsv',
#     resources:
#         mem_mb = lambda wildcards, attempt: 1024 * attempt
#     run:
#         import pandas as pd
#         import json

#         names_explained = [
#             (0, 'sw_score', 'Smith-Waterman score of the match'),
#             (1, 'pct_subs', 'Pct. substitutions in matching region compared to the consensus'),
#             (2, 'pct_del', 'Pct. of bases opposite a gap in the query sequence (deleted bp)'),
#             (3, 'pct_ins', 'Pct. of bases opposite a gap in the repeat consensus (inserted bp)'),
#             (4, 'query_name', 'Query sequence name'),
#             (5, 'query_start', 'Starting position of match in query sequence'),
#             (6, 'query_end', 'Ending position of match in query sequence'),
#             (7, 'query_after_bp', 'bp in query sequence past the ending position of match'),
#             (8, 'is_complement_match', 'True: match is with the Complement of the consensus sequence in the database'),
#             (9, 'repeat_name', 'Name of the matching interspersed repeat'),
#             (10, 'repeat_class', 'Class of the repeat'),
#             (11, 'repeat_start', 'Starting position of match in database sequence (using top-strand numbering)'),
#             (12, 'repeat_end', 'Ending position of match in database sequence'),
#             (13, 'repeat_before_bp', 'bp in (complement of) the repeat consensus sequence prior to beginning of the match'),
#             (14, 'match_ID', 'Consecutive ID of match'),
#             (15, 'is_partly_included', 'True: there is a higher-scoring match whose domain partly (<80%) includes the domain of this match')
#         ]

#         df = pd.read_csv(
#             input.table,
#             header=None,
#             index_col=False,
#             delimiter=r"\s+",
#             skip_blank_lines=True,
#             skiprows=[0,1],
#             comment='#',
#             names=[n[1] for n in names_explained]
#         )

#         name_map = json.loads(open(input.rename, 'r').read())

#         for c in ['query_after_bp', 'repeat_start', 'repeat_before_bp']:
#             df[c] = df[c].apply(lambda x: int(x.strip('()')))
#             df[c] = df[c].astype(int)

#         df['is_complement_match'] = df['is_complement_match'].apply(lambda x: True if x == 'C' else False)
#         df['is_complement_match'] = df['is_complement_match'].astype(bool)
#         df['is_partly_included'] = df['is_partly_included'].apply(lambda x: True if x == '*' else False)
#         df['is_partly_included'] = df['is_partly_included'].astype(bool)
#         df['query_name'].replace(name_map, inplace=True)
#         df.sort_values(['query_name', 'query_start', 'sw_score'], ascending=[True, True, False], inplace=True)

#         with open(output.tsv, 'w') as table:
#             for pos, head, comment in names_explained:
#                 _ = table.write(f'##.{pos} - {head} - {comment}\n')
#             df.to_csv(table, sep='\t', header=True, index=False)
#     # END OF RUN BLOCK
