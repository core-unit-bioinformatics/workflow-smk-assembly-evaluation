

def get_genemodel(gene_model):

    matched_model = config["genemodels"][gene_model]
    model_path = DIR_GLOBAL_REF.joinpath(matched_model)
    return model_path


def get_reference_genome(ref_genome):

    matched_ref = config["refgenomes"][refgenome]["any"]
    ref_path = DIR_GLOBAL_REF.joinpath(matched_ref)
    return ref_path
