#!/bin/bash

module load singularity
eval "$(conda shell.bash hook)"
conda activate snakemake

cd ..
snakemake --snakefile Snakefile --configfile test/config/simulation.yml --use-singularity --cores 1
