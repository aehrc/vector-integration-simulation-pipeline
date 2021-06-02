# Viral integration simulation

Scripts and snakemake workflow for simulating integration of a virus/vector into a host.

## Overview

In order to simulate data, the snakemake will:
1. Parses the config file to generate all combinations of the specified paramters.  Each combination is a 'condition', and each condition has one or more 'replicates'.  All replicates and conditions have a different random seed.
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

Specify the output directory with the key `out_directory`.  Output files can be found in this directory, under the dataset name.  The path should be either absolute, or relative to the snakefile.

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

## Outputs

The main outputs of the pipeline are:
1. Fasta file containing host sequence with integrated viral sequences (and episomes, if appropriate)
2. Paired-end reads from `art_illumina` in fastq and sam format
3. Text file containing location and properties of each integration, in a tab-separated format

The text file with the properties of each integration ('int-info') has the following columns:

 - id: a unique number for each integration
 - chr: the name of the host reference chromosome/contig in which integration occurred
 - hPos: position in the original host reference chromosome/contig at which integration occurred
 - leftStart, leftStop: coordinates of the ambiguous bases (containing a gap or overlap if there is one) on the left side of the integration in the output fasta file.  This will only be the same as hPos for the first integration on each chromsome/contig
 - rightStart, rightStop: same as leftStart, leftStop but for the right side of the integration
 - hDeleted: number of bases deleted from the host chromosome/contig at the integration site.  Deleted bases occur after the right side of the integration
 - hDelted_input_fasta: probably just ignore this column - it should be the same as hDeleted
 - virus: name of viral reference from which integrated sequence was taken
 - vBreakpoints: a list of the parts of the virus which were integrated.  For example [0, 1611] means that bases 0 through 1611 from the viral reference were integrated.  [1400, 1531];[1335, 1400] means that there was a rearrangement: bases 1400 through 1531 and bases 1335 through 1400 were inserted
 - vOris: the orientation of each part of the virus listed in vBreakpoints, + for sense and - for antisense
 - juncTypes: the type of the left and right junctions.  'clean' means there's nothing between the host and viral sequence, 'gap' means random bases were inserted at the jucntion, 'overlap' means that there was homology between the host and vector at the junction.  This is a comma-separated list with two elements - the first is the left junction, and the second is the right junction.
 - juncBases: the sequences at the left and right junctions.  This is a comma separated list with two elements - the first is for the left junction and the second is for the right
 - juncLengths: the number of bases involved in the left and right junctions.This is a comma separated list with two elements - the first is for the left junction and the second is for the right
 - whole: True if the whole virus was integrated, False otherwise
 - rearrangement: True if the integration involved a rearrangement (integration of two or more non-contiguous pieces of the viral reference), False otherwise
 - deletion: True if the part of the virus integrated harbors a deletion, false otherwise
 - n\_swaps: If there was a rearrangement, the number of times that pieces of the integrated part of the virus were swapped
 - n\_delete: If there was a deletion, the number of pieces deleted from the integrated part of the virus

After annotation with simulated reads, the following columns are added:

 - left_chimeric: Chimeric reads that span the left junction
 - right\_chimeric: Chimeric reads that span the left junction
 - left\_discord: Read pairs which straddle the left junction, so there's one read aligned to host and one to virus
 - right\_discord: Read pairs which straddle the right junction, so there's one read aligned to virus and one to host
 - multiple\_discord: Read pairs which span multiple integration junctions, with one read in host and the other in virus/vector.  These should be called as integrations, but it's not clear what the coordinates of the integration should be.
 - fake\_discord: read pairs where both reads map to either host or vector, but span more than one integration. In other words, read pairs where one both reads map to either host or vector, but are interrupted by a sequence of the other type. For example, a pair where one read maps to the virus of one integration, and the other to the virus of the next integration. These should not be called as integrations, but will appear to have a longer template length than expected when mapped
 - discordant\_and\_chimeric: Read pairs which are both discordant (because one read maps to host, and the other to vector, although one of the reads may not be mapped in its entirety), but also one member of the pair is chimeric
 
The text file with the properties of each episomal sequence has the following columns:

- id: a unique number for each episome
- virus: name of viral reference from which episomal sequence was taken
- start, stop: coordinates of viral reference from which episomal sequence was taken
- pieces: a list of the parts of the virus which are episomal.  For example [0, 1611] means that bases 0 through 1611 from the viral reference constitute the episome.  [1400, 1531];[1335, 1400] means that there was a rearrangement: the episome consists of bases 1400 through 1531 and bases 1335 through 1400 
- oris:  the orientation of each part of the virus listed in pieces column, + for sense and - for antisense
- is\_whole: True if the whole virus is episomal, False otherwise
- is\_rearrange: True if the episome is rearranged (episome of two or more non-contiguous pieces of the viral reference), False otherwise
- is\_deletion: True if the episome harbors a deletion, false otherwise
- n\_swaps: If there was a rearrangement, the number of times that pieces of the episome were randomly swapped
- n\_delete: If there was a deletion, the number of pieces deleted from the episome

