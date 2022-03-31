#!/bin/bash
set -e

cd /opt/simvi
eval "$(/bin/micromamba shell hook -s bash)"
 micromamba activate
/opt/conda/bin/snakemake "$@"
