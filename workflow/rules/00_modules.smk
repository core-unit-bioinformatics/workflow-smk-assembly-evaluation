"""
Use this module to list all includes
required for your pipeline - do not
add your pipeline-specific modules
to "commons/00_commons.smk"
"""

include: "00-prepare/settings.smk"
include: "00-prepare/sample_table.smk"
include: "00-prepare/source_setup.smk"

include: "05-preprocess/merge_tag.smk"
include: "05-preprocess/ncbi_fcs.smk"
include: "05-preprocess/remove_contam.smk"

include: "10-asm-align/pyutils.smk"
include: "10-asm-align/align.smk"

include: "20-read-align/align.smk"
include: "20-read-align/filter.smk"
include: "20-read-align/merge.smk"

include: "50-postprocess/asm_chrom_assign.smk"
include: "50-postprocess/asm_ctg_refcov.smk"
include: "50-postprocess/asm_ctg_readcov.smk"

include: "80-statistics/assemblies.smk"
