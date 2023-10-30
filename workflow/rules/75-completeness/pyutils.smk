

def get_gene_model(gene_model):

    matched_model = config["gene_models"][gene_model]
    model_path = DIR_GLOBAL_REF.joinpath(matched_model)
    return model_path


def get_gene_annotation(gene_annotation):

    matched_annotation = config["gene_annotations"][gene_annotation]
    model_path = DIR_GLOBAL_REF.joinpath(matched_annotation)
    return model_path


def get_reference_genome(refgenome):

    matched_ref = config["refgenomes"][refgenome]["any"]
    ref_path = DIR_GLOBAL_REF.joinpath(matched_ref)
    return ref_path
