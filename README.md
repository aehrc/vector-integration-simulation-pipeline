# Viral integration simulation

A [`snakemake`](https://snakemake.readthedocs.io/en/stable/) workflow for simulating integration of a virus/vector into a host.

## Quickstart

The repo includes a host and viral reference for testing your installation. To run this test, you will need to have one of [`conda`](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html), [`docker`](https://docs.docker.com/get-docker/) or [`singularity`](https://sylabs.io/guides/3.5/user-guide/quick_start.html) installed.

### Conda

If you have [`conda`](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) installed, you can run the pipeline with `snakemake`, which automatically downloads dependencies using `conda`.

```
# clone git repo
git clone https://github.com/aehrc/vector-integration-simulation-pipeline.git
cd vector-integration-simulation-pipeline

# create conda environment called 'snakemake' containing snakemake
conda create -n snakemake -c conda-forge -c bioconda snakemake -y
conda activate snakemake

# run with test data
snakemake --configfile test/config/simulation.yml --jobs 1 --use-conda --conda-frontend conda
```

You can find out more about snakemake options in the [snakemake documentation](https://snakemake.readthedocs.io/en/stable/executing/cli.html). 

To run with your own host and viral references, you will need to specify these in the config file (see below).

### Container (`docker` or `singularity`)

Another way to run is using the `docker` container, which contains all dependencies, with either `docker` or `singularity`.  To get the results from inside the container, you will need to [bind-mount](https://docs.docker.com/storage/bind-mounts/) the directory where you would like the results to be written when you run the container.  For example:

```
docker run \
--rm -it \
-v ${PWD}/test_results:/opt/simvi/test/out szsctt/simvi:latest \
snakemake --jobs 1 --configfile test/config/simulation.yml
```

The results should appear in a directory created in the current working directory called `test_results`.  Compare these results to those in the directory `example_results`.

To run with your own host and viral references, you will need to specify these in the config file (see below), and bind-mount the config file and the references into the container.

It's also possible to run the container with `singularity` instead of `docker`, which requires slightly different [bind-mounts](https://sylabs.io/guides/3.5/user-guide/bind_paths_and_mounts.html):

```
singularity pull simvi.sif docker://szsctt/simvi:latest

mkdir -p test_results/.snakemake
singularity exec \
-B ${PWD}/test_results:/opt/simvi/test/out \
-B ${PWD}/test_results/.snakemake:/opt/simvi/.snakemake \
simvi.sif \
/opt/simvi/run_singularity.sh --configfile /opt/simvi/test/config/simulation.yml --jobs 1
```

You can also give additional arugments which will be passed to `snakemake`.

Note that since singularity containers are read-only, we need to also bind-mount a folder for snakemake to write to (in this case, `/opt/simvi/.snakemake` inside the container to `${PWD}/test_results/.snakemake` outside the container).

## Overview

In order to simulate data, `snakemake` will:
1. Parses the config file to generate all combinations of the specified paramters.  Each combination of parameters is a 'condition', and each condition has one or more 'replicates'.  All replicates and conditions have a different random seed.
2. Simulates integration using the `python3` script `scripts/insert_virus.py`.  This script outputs a `fasta` file consisting of the host with viral sequences inserted, and two files describing the location and properties of the integrated virus and episomes
3. `art_illumina` is used to simulate paried-end reads based on the `fasta` file from the previous step
4. A script `scripts/annotate_reads.py` annotates the reads crossing host/viral junctions for each integration

For each dataset, the reads can be found in the output directory under `sim_reads`, and the information about the integrations can be found under `sim_ints.`

## Dependencies

As outlined in the [*QuickStart*](https://github.com/aehrc/vector-integration-simulation-pipeline/edit/master/README.md#quickstart) section above, running the pipeline requires either `snakemake` and `conda`/`mamba`, or `docker`/`singularity`.  If using `conda`, other dependencies are automatically downloaded by `snakemake` (using the `--use-conda`) argument using the environment files in the `envs` directory.  If using the container with `docker`/`singularity`, dependences are already installed inside the container.


## Config file

Specify all inputs and options in a config file, which is provided to `snakemake`.  An example config file can be found at `test/config/simulation.yml`.

The config file is organised into datasets - in the example above, there is one dataset `test`.  A config file should have one or more datasets. The results for each dataset will be output in a directory with the same name as the dataset (`test` in the case of the supplied config).

If you use `docker` or `singularity` to run the pipeline with your own references, make sure to [bind-mount](https://docs.docker.com/storage/bind-mounts/) your data and output directories into the container, and that any paths in your config file are correct inside the container.

For each dataset, the following parameters should be specified in the config file.

#### Output directory

Specify the output directory with the key `out_directory`.  Output files can be found in this directory, under the dataset name.  The path should be either absolute, or relative to the directory in which you run the workflow (usually the directory which you cloned the repo into).  If you are using the container version, the `Snakefile` is located at `/opt/simvi/Snakefile`.

#### Host, viral references

Use the keys `hosts` and `viruses` to specify a dict of host and viruses, respectivley, in which the keys in each dict are a name for that reference, and the value is the path (absolute or relative to the snakefile) to that reference in `fasta` format.  

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

The text file with the properties of each integration ('*int-info.tsv') has the following columns:

 - `id`: a unique number for each integration
 - `chr`: the name of the host reference chromosome/contig in which integration occurred
 - `hPos`: position in the original host reference chromosome/contig at which integration occurred
 - `leftStart`, `leftStop`: coordinates of the ambiguous bases (containing a gap or overlap if there is one) on the left side of the integration in the output fasta file.  This will only be the same as `hPos` for the first integration on each chromsome/contig
 - `rightStart`, `rightStop`: same as `leftStart`, `leftStop` but for the right side of the integration
 - `hDeleted`: number of bases deleted from the host chromosome/contig at the integration site.  Deleted bases occur after the right side of the integration
 - `hDelted_input_fasta`: probably just ignore this column - it should be the same as `hDeleted`
 - `virus`: name of viral reference from which integrated sequence was taken
 - `vBreakpoints`: a list of the parts of the virus which were integrated.  For example [0, 1611] means that bases 0 through 1611 from the viral reference were integrated.  [1400, 1531];[1335, 1400] means that there was a rearrangement: bases 1400 through 1531 and bases 1335 through 1400 were inserted
 - `vOris`: the orientation of each part of the virus listed in vBreakpoints, + for sense and - for antisense
 - `juncTypes`: the type of the left and right junctions.  'clean' means there's nothing between the host and viral sequence, 'gap' means random bases were inserted at the jucntion, 'overlap' means that there was homology between the host and vector at the junction.  This is a comma-separated list with two elements - the first is the left junction, and the second is the right junction.
 - `juncBases`: the sequences at the left and right junctions.  This is a comma separated list with two elements - the first is for the left junction and the second is for the right
 - `juncLengths`: the number of bases involved in the left and right junctions.This is a comma separated list with two elements - the first is for the left junction and the second is for the right
 - `whole`: True if the whole virus was integrated, False otherwise
 - `rearrangement`: True if the integration involved a rearrangement (integration of two or more non-contiguous pieces of the viral reference), False otherwise
 - `deletion`: True if the part of the virus integrated harbors a deletion, false otherwise
 - `n\_swaps`: If there was a rearrangement, the number of times that pieces of the integrated part of the virus were swapped
 - `n\_delete`: If there was a deletion, the number of pieces deleted from the integrated part of the virus

After annotation with simulated reads, the following columns are added  ('*int-info.annotated.tsv'):

 - `left_chimeric`: Read IDs of the chimeric reads that span the left junction
 - `right\_chimeric`: Read IDs of the  chimeric reads that span the left junction
 - `left\_discord`: Read pairs which straddle the left junction, so there's one read aligned to host and one to virus
 - `right\_discord`: Read pairs which straddle the right junction, so there's one read aligned to virus and one to host
 - `multiple\_discord`: Read pairs which span multiple integration junctions, with one read in host and the other in virus/vector.  These should be called as integrations, but it's not clear what the coordinates of the integration should be.
 - `fake\_discord`: read pairs where both reads map to either host or vector, but span more than one integration. In other words, read pairs where one both reads map to either host or vector, but are interrupted by a sequence of the other type. For example, a pair where one read maps to the virus of one integration, and the other to the virus of the next integration. These should not be called as integrations, but will appear to have a longer template length than expected when mapped
 - `discordant\_and\_chimeric`: Read pairs which are both discordant (because one read maps to host, and the other to vector, although one of the reads may not be mapped in its entirety), but also one member of the pair is chimeric
 
The text file with the properties of each episomal sequence ('*epi-info.tsv') has the following columns:

- `id`: a unique number for each episome
- `virus`: name of viral reference from which episomal sequence was taken
- `start`, `stop`: coordinates of viral reference from which episomal sequence was taken
- `pieces`: a list of the parts of the virus which are episomal.  For example [0, 1611] means that bases 0 through 1611 from the viral reference constitute the episome.  [1400, 1531];[1335, 1400] means that there was a rearrangement: the episome consists of bases 1400 through 1531 and bases 1335 through 1400 
- `oris`:  the orientation of each part of the virus listed in pieces column, + for sense and - for antisense
- `is\_whole`: True if the whole virus is episomal, False otherwise
- `is\_rearrange`: True if the episome is rearranged (episome of two or more non-contiguous pieces of the viral reference), False otherwise
- `is\_deletion`: True if the episome harbors a deletion, false otherwise
- `n\_swaps`: If there was a rearrangement, the number of times that pieces of the episome were randomly swapped
- `n\_delete`: If there was a deletion, the number of pieces deleted from the episome

