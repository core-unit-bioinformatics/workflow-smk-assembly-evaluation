
PATH_ID_LENGTH = 8

DATA_ROOT = config.get("data_root", "/")
DATA_ROOT = pathlib.Path(DATA_ROOT).resolve(strict=True)

WILDCARDS_REF_GENOMES = ["hg38", "t2tv2"]

# Tool parameter

## mosdepth

### window size for assembly-to-reference
### alignments; to evaluate assembly contig
### coverage in reference genome

_tmp = config.get("mosdepth", dict())
MOSDEPTH_ASSM_REF_COV_WINDOW_SIZE = int(_tmp.get("assm_ref_cov_window_size", 10000))
_tmp = None
