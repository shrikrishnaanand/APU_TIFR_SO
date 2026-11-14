#!/bin/bash
#SBATCH --gres=gpu:1
#SBATCH --partition=gpu
#SBATCH --time=24:00:00
#SBATCH --job-name=et_e0_5
#SBATCH --output=et_e0_5_%j.out
#SBATCH --error=et_e0_5_%j.out

# Usage (from inside et_e0/):  sbatch submit_et_e0_5.sh
# The notebook uses relative paths (NLO_e0/, NNLO_e0/, ../x_e_files/),
# so everything runs from the et_e0 directory.

cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
echo "== workdir: $(pwd)"

# Activate environment (uenv may print harmless libtinfo/GLIBCXX warnings;
# don't let them abort the job)
uenv verbose cuda-11.4.4 cudnn-11.x-8.8.0 || true
uenv miniconda3-py39 || true
conda activate sage || { echo "!! conda activate sage failed"; exit 1; }

# uenv puts /opt/miniconda39/lib on LD_LIBRARY_PATH, which shadows the system
# libtinfo/libstdc++ and breaks system tools (bash, dig, ...). The conda env
# finds its own libs via RPATH, so drop that entry (CUDA paths are kept).
export LD_LIBRARY_PATH=$(echo "${LD_LIBRARY_PATH:-}" | tr ':' '\n' | grep -v '^/opt/miniconda39' | paste -sd: -)

echo "== python:  $(which python)"
echo "== sage:    $(which sage)"
echo "== jupyter: $(which jupyter)"

# From here on, stop at the first real error
set -eo pipefail

NB=et_e0_5.ipynb
BASE=${NB%.ipynb}
mkdir -p NNLO_e0

echo "== 1) nbconvert $NB -> $BASE.sage"
jupyter nbconvert --to script "$NB" --stdout > "$BASE.sage"

echo "== 2) sage --preparse $BASE.sage -> $BASE.sage.py"
sage --preparse "$BASE.sage"

echo "== 3) running $BASE.sage.py"
export PYTHONUNBUFFERED=1
sage "$BASE.sage.py"
echo "== done"
