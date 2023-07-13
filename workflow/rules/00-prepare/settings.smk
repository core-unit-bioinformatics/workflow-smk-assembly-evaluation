
PATH_ID_LENGTH = 8

DATA_ROOT = config.get("data_root", "/")
DATA_ROOT = pathlib.Path(DATA_ROOT).resolve(strict=True)

WILDCARDS_REF_GENOMES = ["hg38", "t2tv2"]
