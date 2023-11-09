
localrules: estimate_asm_unit_karyotype
rule estimate_asm_unit_karyotype:
    """
    NB: if the output path is changed, may need to
    adapt code in 75-completeness::pyutils !
    """
    input:
        ref_assign = expand(
            rules.estimate_chromosome_assignment.output.tsv_target,
            asm_unit=ASSEMBLY_UNITS_SEX_SPECIFIC,
            allow_missing=True
        )
    output:
        karyo_est = DIR_RES.joinpath(
            "reports", "ref_chrom_assign",
            "{sample}.{refgenome}.karyo-est.tsv"
        )
    run:
        import pathlib as pl
        import pandas as pd

        force_assign_any = False
        sex_chromosomes = config.get("sex_chromosomes", dict())
        if FORCE_ANNOTATED_SAMPLE_SEX:
            pass
        elif not sex_chromosomes:
            logerr(
                (
                    "No config entry 'sex_chromosomes' found - "
                    "assigning karyotype 'any' to all assembly units."
                ))
            force_assign_any = True
        else:
            if not all(label in sex_chromosomes for label in ["male", "female"]):
                logerr(f"Sex chromosomes must be specified for both sexes (male/female): {sex_chromosomes}")
                raise ValueError(f"Config entry 'sex_chromosomes' is malformed: {sex_chromosomes}")

        out_records = []
        for tsv_file in input.ref_assign:
            sample, asm_unit, refgenome, _, tsv_ext = pl.Path(tsv_file).name.rsplit(".", 4)
            assert tsv_ext == "tsv"
            if FORCE_ANNOTATED_SAMPLE_SEX:
                karyotype = SAMPLE_INFOS[sample]["sex"]
                out_records.append((sample, asm_unit, refgenome, karyotype, -1, -1, -1, -1))
                continue
            if force_assign_any:
                karyotype = "any"
                out_records.append((sample, asm_unit, refgenome, karyotype, -1, -1, -1, -1))
                continue
            df = pd.read_csv(tsv_file, sep="\t", header=0, comment="#")
            female_score, female_frag = compute_karyotype_score(df, sex_chromosomes["female"])
            male_score, male_frag = compute_karyotype_score(df, sex_chromosomes["male"])
            if female_score > male_score:
                karyotype = "female"
            elif male_score > female_score:
                karyotype = "male"
            elif male_score == female_score and male_score > 0:
                logerr(
                    (
                        f"Karyotype scores identical for sample / assembly unit: {sample} / {asm_unit}\n"
                        f"Female: {female_score} /// Male: {male_score} --- assigning 'any'."
                    ))
                karyotype = "any"
            else:
                # both scores zero
                logerr(f"All karyotype scores are zero - check the alignments? Assigning karyotype 'any'.")
                karyotype = "any"
            out_records.append(
                (sample, asm_unit, ref, karyotype,
                female_score, female_frag, male_score, male_frag)
            )

        out_records = pd.DataFrame.from_records(out_records,
            columns=[
                "sample", "asm_unit", "refgenome", "karyotype",
                "female_score", "female_fragmentation",
                "male_score", "male_fragmentation"
            ]
        )
        out_records.sort_values("asm_unit", inplace=True)
        out_records.to_csv(output.karyo_est, sep="\t", header=True, index=False)
    # END OF RUN BLOCK


localrules: build_assembly_karyotype_summary
rule build_assembly_karyotype_summary:
    input:
        all_est = expand(
            rules.estimate_asm_unit_karyotype.output.karyo_est,
            sample=SAMPLES,
            refgenome=COMPLETE_REF_GENOME
        )
    output:
        tsv = DIR_RES.joinpath(
            "reports", "ref_chrom_assign",
            "all-samples.karyo-est.tsv"
        )
    run:
        import pandas as pd

        merged = []
        for tsv_file in input.all_est:
            merged.append(pd.read_csv(tsv_file, sep="\t", header=0))
        merged = pd.concat(merged, axis=0, ignore_index=False)
        merged.sort_values(["sample", "asm_unit"], inplace=True)
        merged.to_csv(output.tsv, sep="\t", header=True, index=False)
