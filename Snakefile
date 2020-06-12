from os import path


#####################################################
################### check inputs ####################
#####################################################

# check required inputs: require host and vector fasta
if 'hosts' not in config or 'viruses' not in config:
	raise ValueError("Please specify 'hosts' and 'viruses' in config file")
	
if len(config['hosts']) == 0 or len(config['viruses']) == 0:
	raise ValueError("length of host or viruses dictionary is zero. Please specify at least one file to use")


def check_and_set_default(field, default, name=""):
	if field not in config:
		config[field] = default
		print(f"{name} ({field}) not set, using {default}")

# check for output directory and normalise this path
check_and_set_default('out_directory', '.', "output directory")
config['out_directory'] = path.normpath(config['out_directory'])

# check for number of replicates (default is 1)
check_and_set_default('replicates', 1, "number of replicates")

# check parameters for integration simulation
check_and_set_default('initial_seed', 42, "random seed for first replicate")
check_and_set_default('num_ints', 5, "number of integrations")
	
# check art-illumina settings
check_and_set_default('read_len', 150, "read length of simulated reads")
check_and_set_default('fcov', 10, "fold coverage of simulated reads")
check_and_set_default('frag_len', 500, "mean fragment length for simulated reads")
check_and_set_default('frag_std', 30, "fragment length standard deviation for simulated reads")
check_and_set_default('seq_sys', "HS25", "sequencing system")
check_and_set_default('art_seed', 42, "random seed")


#####################################################
################### target files ####################
#####################################################

rule all:
	input: 
		expand("{outdir}/reads/sim_ints.{host}.{virus}.{replicate}.sam", outdir = config['out_directory'], host = list(config['hosts'].keys()), virus = list(config['viruses'].keys()), replicate = [int(i) for i in range(config['replicates'])])
		


#####################################################
############### simulate integrations ###############
#####################################################


rule simulate_integrations:
	input:
		host = lambda wildcards: config['hosts'][wildcards.host],
		virus =  lambda wildcards: config['viruses'][wildcards.virus]
	output:
		sim_fasta = "{outdir}/integrations/sim_ints.{host}.{virus}.{replicate}.fa",
		sim_info = "{outdir}/integrations/sim_info.{host}.{virus}.{replicate}.tsv"
	params:
		int_num = config['num_ints'],
		seed = lambda wildcards: int(config['initial_seed']) + int(wildcards.replicate)
	shell:
		"""
		python3 insert_virus_simple.py \
		--host {input.host} \
		--virus {input.virus} \
		--ints {output.sim_fasta} \
		--info {output.sim_info} \
		--int_num {params.int_num} \
		--seed {params.seed}
		"""

#####################################################
################### simulate reads ##################
#####################################################

rule art:
	input:
		sim_fasta = rules.simulate_integrations.output.sim_fasta,
	output:
		sam = "{outdir}/reads/sim_ints.{host}.{virus}.{replicate}.sam"
	params:
		seq_system = f"-ss {config['seq_sys']}",
		read_length = f"-l {config['read_len']}",
		coverage = f"-f {config['fcov']}",
		frag_length = f"-m {config['frag_len']}",
		frag_std = f"-s {config['frag_std']}",
		seed = lambda wildcards: f"-rs {int(config['art_seed']) + int(wildcards.replicate)}",
		paried = "-p",
		sam = "-sam",
		input = lambda wildcards, input: f"-i {input.sim_fasta}",
		output = lambda wildcards, output: f"-o {path.splitext(output.sam)[0]}"
	shell:
		"""
		art_illumina {params}
		"""
		
		
			
		