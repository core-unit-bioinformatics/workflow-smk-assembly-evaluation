
PATH_ID_LENGTH = 8

DATA_ROOT = config.get("data_root", "/")
DATA_ROOT = pathlib.Path(DATA_ROOT).resolve(strict=True)

WILDCARDS_REF_GENOMES = ["hg38", "t2tv2"]

# STATISTICS PARAMETERS

SEQUENCE_LENGTH_THRESHOLDS_ASSEMBLY = {
    "default": [10000, 100000, 1000000, 10000000]
}
for asm_unit, thresholds in config.get("sequence_length_thresholds_assembly", dict()).items():
    if asm_unit.startswith("setting"):
        continue
    elif asm_unit.startswith("contaminants"):
        lookup_name = asm_unit
    elif asm_unit.startswith("asm"):
        lookup_name = asm_unit
    else:
        lookup_name = f"asm_{asm_unit}"
    SEQUENCE_LENGTH_THRESHOLDS_ASSEMBLY[lookup_name] = thresholds
ASSEMBLY_UNITS_RELEVANT_DEFAULT = sorted(k for k in SEQUENCE_LENGTH_THRESHOLDS_ASSEMBLY.keys() if k != "default")

# TOOL PARAMETERS

## mosdepth

### mosdepth global on/off switch
RUN_MOSDEPTH = config.get("run_mosdepth", False)

### window size for assembly-to-reference
### alignments; to evaluate assembly contig
### coverage in reference genome
RUN_MOSDEPTH_ASSM_REF_COV = RUN_MOSDEPTH and config.get("run_mosdepth_assm_ref_cov", False)
MOSDEPTH_PARAMS = config.get("mosdepth", dict())
MOSDEPTH_ASSM_REF_COV_WINDOW_SIZE = int(MOSDEPTH_PARAMS.get("assm_ref_cov_window_size", 10000))


## NCBI FCS toolkit

RUN_NCBI_FCS_ADAPTOR = config.get("run_ncbi_fcs_adaptor", False)
RUN_NCBI_FCS_GX = config.get("run_ncbi_fcs_gx", False)

NCBI_FCS_PARAMS = config.get("ncbi_fcs", dict())
NCBI_FCS_ROOT_PATH = pathlib.Path(NCBI_FCS_PARAMS.get("root_path", "/"))
NCBI_FCS_ROOT_PATH = NCBI_FCS_ROOT_PATH.resolve(strict=True)

NCBI_FCS_ADAPTOR_SIF = pathlib.Path("CONFIG.NOT-SET.NCBI_FCS_ADAPTOR_SIF")
NCBI_FCS_ADAPTOR_SCRIPT = pathlib.Path("CONFIG.NOT-SET.NCBI_FCS_ADAPTOR_SCRIPT")
NCBI_FCS_ADAPTOR_TAXONOMY = ""
if RUN_NCBI_FCS_ADAPTOR:
    NCBI_FCS_ADAPTOR_SIF = NCBI_FCS_ROOT_PATH.joinpath(
        NCBI_FCS_PARAMS["adaptor_sif"]
    )
    assert NCBI_FCS_ADAPTOR_SIF.is_file()
    NCBI_FCS_ADAPTOR_SCRIPT = NCBI_FCS_ROOT_PATH.joinpath(
        NCBI_FCS_PARAMS["adaptor_script"]
    )
    assert NCBI_FCS_ADAPTOR_SCRIPT.is_file()
    NCBI_FCS_ADAPTOR_TAXONOMY = NCBI_FCS_PARAMS["adaptor_taxonomy"]


NCBI_FCS_GX_DB_PATH = pathlib.Path("CONFIG.NOT-SET.NCBI_FCS_GX_DB_PATH")
NCBI_FCS_GX_DB_NAME = NCBI_FCS_PARAMS.get("gx_db_name", "CONFIG.NOT-SET.NCBI_FCS_GX_DB_NAME")
NCBI_FCS_GX_SIF = pathlib.Path("CONFIG.NOT-SET.NCBI_FCS_GX_SIF")
NCBI_FCS_GX_SCRIPT = pathlib.Path("CONFIG.NOT-SET.NCBI_FCS_GX_SCRIPT")
NCBI_FCS_GX_TAX_ID = 0
if RUN_NCBI_FCS_GX:
    NCBI_FCS_GX_DB_PATH = NCBI_FCS_ROOT_PATH.joinpath(
        NCBI_FCS_PARAMS["gx_db_path"]
    )
    assert NCBI_FCS_GX_DB_PATH.is_dir()
    assert not NCBI_FCS_GX_DB_NAME.startswith("CONFIG.NOT-SET")
    NCBI_FCS_GX_SIF = NCBI_FCS_ROOT_PATH.joinpath(
        NCBI_FCS_PARAMS["gx_sif"]
    )
    assert NCBI_FCS_GX_SIF.is_file()
    NCBI_FCS_GX_SCRIPT = NCBI_FCS_ROOT_PATH.joinpath(
        NCBI_FCS_PARAMS["gx_script"]
    )
    assert NCBI_FCS_GX_SCRIPT.is_file()
    NCBI_FCS_GX_TAX_ID = NCBI_FCS_PARAMS["gx_tax_id"]
