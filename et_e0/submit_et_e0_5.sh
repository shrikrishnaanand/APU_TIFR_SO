#!/bin/bash
#SBATCH --gres=gpu:1
#SBATCH --partition=gpu
#SBATCH --time=24:00:00
#SBATCH --job-name=et_e0_5
#SBATCH --output=et_e0_5_%j.out

# Usage (from inside et_e0/):  sbatch submit_et_e0_5.sh
# The notebook uses relative paths (NLO_e0/, NNLO_e0/, ../x_e_files/),
# so everything runs from the et_e0 directory.

set -eo pipefail
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"

# Activate environment
uenv verbose cuda-11.4.4 cudnn-11.x-8.8.0
uenv miniconda3-py39
conda activate sage

NB=et_e0_5.ipynb
BASE=${NB%.ipynb}
mkdir -p NNLO_e0

# 1) Notebook -> script, renamed to .sage
jupyter nbconvert --to script "$NB" --stdout > "$BASE.sage"

# 2) .sage -> .sage.py  (Sage preparser: handles ^, 1/2 rationals, var(), etc.)
sage --preparse "$BASE.sage"

# 3) Run it (unbuffered so progress shows up in the .out file)
export PYTHONUNBUFFERED=1
sage "$BASE.sage.py"
