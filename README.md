# Viral integraiton simulation

Scripts and snakemake workflow for simulating integration of a virus/vector into a host.

## Overview

In order to simulate data, the snakemake will:
1. Parses the config file to generate all combinations of the specified paramters.  Each combionation is a 'condition', and each condition has one or more 'replicates'
2. Simulates integration using the `python3` script `scripts/insert_virus.py`.  This script outputs a `fasta` file consisting of the host with viral sequences inserted, and two files describing the location and properties of the integrated virus and episomes
3. `art_illumina` is used to simulate paried-end reads based on the `fasta` file from the previous step
4. A script `scripts/annotate_reads.py` annotates the reads crossing host/viral junctions for each integration

For each dataset, the reads can be found in the output directory under `sim_reads`, and the information about the integrations can be found under `sim_ints.`


## Dependencies

Execution requires `snakemake` and either `singularity` (recommended) or `conda` (envs may be out of date) to supply other dependencies.  If running with the `--use-singularity` flag, you may need to bind-mount the directories containing data and ouputs via the `--singularity-args` flag (eg `--singularity-args '-B /path/to/data -B /path/to/output'`).

## Test script

To run a demo, use `test/runMe.sh`.  Note that this script requires the user to have a conda environment called `snakemake`, with `snakemake` installed, and a modulefile called `singularity`.  Adjust accordingly if your setup is not the same.

## Inputs

Specify all inputs in a config file, which is provided to snakemake.  An example config file is as follows:

```
test:
    out_directory: "test/out/" 
    replicates: 2 

    hosts: 
        "human" : "test/refs/test_human.fa"
    viruses:
        "AAV" :  "test/refs/test_AAV.fa"

    initial_seed: 12345
    seed_increment: 123 
    
    int_num: 
        - 3
    min_sep: 
        - 50  
    epi_num:
        - 0

    p_whole:
        - 0.5
    min_len: 
        - 50 
    p_rearrange:
        - 0.1
    p_delete:
        - 0.1
    lambda_split:
        - 1
    p_overlap:
        - 0.2
    p_gap:
        - 0.2
    lambda_junction:
        - 1
    p_host_deletion:
        - 0.2
    lambda_host_deletion:
        - 20

    read_len: 
        - 150
    fcov: 
        - 10 
    frag_len: 
        - 250
        - 500 
    frag_std: 
        - 30
    seq_sys: 
        - "MSv3"
```

The config file is organised into datasets - in the example above, there is one dataset `test`.  A config file should have one or more datasets.

For each dataset, the following parameters should be specified:

#### Output directory

Specify the output directory with the key `out_directory`.  Output files can be found in this directory, under the dataset name.  The path should be either absolute or relative to the snakefile.

#### Host, viral references

Use the keys `hosts` and `viruses` to specify a dict of host and viruses, respectivley, in which the keys in each dict are a name for that reference, and the value is the path (absolute or relative to the snakefile) to that reference.  

#### Replicates

Each unique combination of simulation parameters (including the host and viral references specified above) is a 'condition'.  Specify the number of replicates to perform for each condition.  Each replicate will have a random seed - specify the seed for the first condition with the key `initial_seed`, and each additional replicate will have a seed incremented by the key `seed_increment`.  This random seed is used for both simuating integrations (python script), and simulating reads (`art_illumina`).

#### Integration properties

The user may specify a number of desired integration properties.  Most properties are specified as either probabilities (start with `p_`) or the mean of a Poisson distribution (`lambda_`).

- `int_num`: The number of integrations
- `min_sep`: The minimum separation (in bp) between each integration
- `epi_num`: The number of extra viral sequences (episomes) added to the output fasta.  Episomes are also subject to being whole/subsequences, rearrangements and deletions
- `p_whole`: Probability that each integration/episome will consist of the whole viral sequence.
- `min_len`: If integration/episome is not whole, its minimum length (in bp)
- `p_rearrange`: Probability that an integration/episome is rearranged
- `p_delete`: Probability that an integration/episome has a deletion
- `lambda_split`: During rearrangmenet/deletion, number of fragments in which to split the viral sequence
- `p_overlap`: Probability of an overlap (common sequence between host and virus) at each junction
- `p_gap`: Probability of a gap (random bases added between host and virus) at each junction
- `lambda_junction`: Mean of Poisson distribution of number of bases involved in each gap/overlap junction
- `p_host_deletion`: Probability of a deletion from the host at the integration site
- `lambda_host_deletion`: Mean of Poisson distribution of number of bases deleted from host at integration site

#### `art_illumina` parameters

The parameters `read_len` (read length), `fcov` (fold coverage), `frag_len` (mean fragment length), `frag_std` (standard deviation of fragment length), and `seq_sys` will be used during read simulation with `art_illumina`.  Further details of these paramters can be found at the [art manpage](https://manpages.debian.org/stretch/art-nextgen-simulation-tools/art_illumina.1.en.html).


