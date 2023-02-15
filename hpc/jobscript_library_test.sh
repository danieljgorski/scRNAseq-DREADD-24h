#!/usr/bin/bash

#PBS -A ECM-Cardiac-Ischemia
#PBS -N scRNAseq-DREADD-24h
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=1:mem=5G
#PBS -m abe
#PBS -M gorskid@hhu.de

# ignore this line, but always keep it
set +eu

WORKDIR=/gpfs/project/gorskid/scRNAseq-DREADD-24h

# change to workdir
cd ${WORKDIR}

# make conda available
module load Miniconda/3

# activate conda env
conda activate r4.2

# ignore this line, but always keep it
set -euo pipefail

# say, execute your R script
# The script needs to be located in the workdir
Rscript hpc/libraries.R