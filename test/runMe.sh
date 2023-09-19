#!/bin/bash

snakemake --snakefile Snakefile --configfile test/config/simulation.yml --scheduler greedy --cores 1
