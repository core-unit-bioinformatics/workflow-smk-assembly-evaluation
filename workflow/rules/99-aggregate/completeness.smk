
ASSEMBLY_COMPLETENESS_OUTPUT = []

if True:  # TODO: make proper switch

    ASSEMBLY_COMPLETENESS_OUTPUT.extend(
        rules.run_all_asm_completeness_asmgene.input.asmgene_stats
    )
