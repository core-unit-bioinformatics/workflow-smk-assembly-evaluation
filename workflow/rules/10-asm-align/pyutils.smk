
def get_reference_assembly(sample, reference):

    sample_sex = SAMPLE_INFOS[sample]["sex"]
    matched_ref = config["refgenomes"][reference][sample_sex]
    ref_path = DIR_GLOBAL_REF.joinpath(matched_ref)
    return ref_path
