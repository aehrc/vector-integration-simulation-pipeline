#!/bin/bash

snakemake --snakefile Snakefile --configfile test/config/simulation.yml --cores 1
