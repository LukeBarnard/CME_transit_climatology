#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --threads-per-core=1
#SBATCH --partition=long

#SBATCH --job-name=cme_transit_climatology
#SBATCH --output=cme_transit_climatology_report.txt

#SBATCH --time=50:00:00

module load anaconda
source activate sir_huxt
python article.py
