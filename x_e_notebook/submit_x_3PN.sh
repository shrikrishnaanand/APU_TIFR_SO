#!/bin/bash
#SBATCH --job-name=pn_2_15
#SBATCH --output=pn_2_15_%j.out
#SBATCH --error=pn_2_15_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --time=24:00:00
#SBATCH --mem=32G

# Load/activate Miniforge
source ~/miniforge3/bin/activate

# Activate Sage environment
conda activate sage

# Run Sage script
sage x_3PN.sage.py
