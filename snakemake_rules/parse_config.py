from os import path
import pandas as pd
import itertools

## default values to fill in if they aren't set
# random seed - for both python script and ART
default_replicates = 3
default_initial_seed = 1234567
default_seed_increment = 5

# genome variation parameters
default_snp_count = 100
default_indel_count = 10
default_cnv_count = 3
default_inversion_count = 1
default_translocation_count = 0

# integration parameters
default_seed_increment = 1000
default_int_num = [5]
default_p_whole = [0.5]
default_p_rearrange = [0.05]
default_p_delete = [0.05]
default_lambda_split = [1.5]
default_p_overlap = [0.3]
default_p_gap = [0.3]
default_lambda_junction = [4]
default_p_host_deletion = [0.2]
default_lambda_host_deletion = [1000]
default_min_sep = [500]
default_min_len = [500]

# read simulation parameters
default_read_len = [150]
default_fcov = [10]
default_frag_len = [500]
default_frag_std = [30]
default_seq_sys = ["HS25"]


def parse_config(config):
	
	##### check required fields for each experiment ####

	# the only requirements are an output directory ('out_directory'), and at least one
	# host and virus fasta file ('hosts' and 'viruses' dicts)
	
	for exp in config:
	
		# check output directory
		check_exists_type_length(config, exp, 'out_directory', str)
		config[exp]['out_directory'] = path.normpath(config[exp]['out_directory'])

		# check host reference(s)
		check_exists_type_length(config, exp, 'hosts', dict)

		# check virus reference(s)
		check_exists_type_length(config, exp, 'viruses', dict)



	#### assign default values if they aren't already assigned ####


	for exp in config:

		# number of replicates
		config = check_and_set_default(config, exp, 'replicates', default_replicates)
		check_required(config, exp, 'replicates')
		check_type(config, exp, 'replicates', int)
		if config[exp]['replicates'] < 1:
			print(f"must have at least one replicate for experiment {exp}: setting number of replicates to 1")
	
		# initial random seed and increment
		config = standard_checks_int(config, exp, 'initial_seed', default_initial_seed)
		config = standard_checks_int(config, exp, 'seed_increment', default_seed_increment)
		
		# genome variation parameters
		config = standard_checks_int(config, exp, 'snp_count', default_snp_count)
		config = standard_checks_int(config, exp, 'indel_count', default_indel_count)		
		config = standard_checks_int(config, exp, 'cnv_count', default_cnv_count)
		config = standard_checks_int(config, exp, 'inversion_count', default_inversion_count)		
		config = standard_checks_int(config, exp, 'translocation_count', default_translocation_count)

		# integration parameters
		config = standard_checks(config, exp, 'min_sep', default_min_sep, list, int) 
		config = standard_checks(config, exp, 'min_len', default_min_len, list, int) 
		config = standard_checks(config, exp, 'int_num', default_int_num, list, int) 
		config = standard_checks(config, exp, 'epi_num', default_int_num, list, int) 
		config = standard_checks(config, exp, 'p_whole', default_p_whole, list, (int, float)) 
		config = standard_checks(config, exp, 'p_rearrange', default_p_rearrange, list, (int,  float)) 
		config = standard_checks(config, exp, 'p_delete', default_p_delete, list, (int, float)) 
		config = standard_checks(config, exp, 'lambda_split', default_lambda_split, list, (int,  float)) 
		config = standard_checks(config, exp, 'p_overlap', default_p_overlap, list, (int,  float)) 
		config = standard_checks(config, exp, 'p_gap', default_p_gap, list, (int,  float)) 
		config = standard_checks(config, exp, 'lambda_junction', default_lambda_junction, list, (int, float)) 
		config = standard_checks(config, exp, 'p_host_deletion', default_p_host_deletion, list, (int, float)) 
		config = standard_checks(config, exp, 'lambda_host_deletion', default_lambda_host_deletion, list, (int, float)) 

		# art parameters
		config = standard_checks(config, exp, 'read_len', default_read_len, list, int) 
		config = standard_checks(config, exp, 'fcov', default_fcov, list, (int, float))
		config = standard_checks(config, exp, 'frag_len', default_frag_len, list, (int, float)) 
		config = standard_checks(config, exp, 'frag_std', default_frag_std, list, (int, float)) 
		config = standard_checks(config, exp, 'seq_sys', default_seq_sys, list, str) 


	##### check for clashes between references for experiments #####
	# reference names and fastas must have a 1-to-1 mapping, and there can't be any 
	# references with the same name or fasta

	# we also need a dict for keeping track of which fasta files correspond to which 
	# reference names

	def add_to_ref_dict(config_refs, ref_dict):
		for name, fasta in config_refs.items():
			if name not in ref_dict:
				ref_dict[name] = fasta
			else:
				if ref_dict[name] != fasta:
					raise ValueError(f"Each name/fasta pair must be unique for all the hosts and references and all experiments.  Found two different fasta files for reference {name}")

	ref_dict = {}
	for exp in config:
		add_to_ref_dict(config[exp]['hosts'], ref_dict)
		add_to_ref_dict(config[exp]['viruses'], ref_dict)


	#### generate experimental conditions and output filename ####

	# within each experiment, we want to run every combination of all the specified parameters
	# each unique combination is a 'condition', 
	# and we want to run it a number of times (the specified number of replicates)

	# give each condition a name, and construct a data frame where each row corresponds to one condition
	# so we can track which parameters we should use for each condition in each experiment


	col_names = [# these three columns should uniquely identify each row
				 'experiment',
				 'condition',
				 'replicate',
				 
				 # parameters generated from product of those specified in config[experiment]
				 # each unique combination in an experiment is one condition 
				 'host_name',
				 'virus_name',
				 'int_num',
				 'epi_num',
				 'min_sep',
				 'min_len',
				 'p_whole',
				 'p_rearrange',
				 'p_delete',
				 'lambda_split',
				 'p_overlap',
				 'p_gap',
				 'lambda_junction',
				 'p_host_deletion',
				 'lambda_host_deletion',
				 'read_len',
				 'fcov',
				 'frag_len',
				 'frag_std',
				 'seq_sys',
				 
				 # columns which depend on the columns above or are constant
				 'out_directory',
				 'host_fasta',
				 'virus_fasta',
				 'random_seed',
				 'snp_count',
				 'indel_count',
				 'cnv_count',
				 'inversion_count',
				 'translocation_count',
				 'sim_fa_filename',
				 'sim_int_info_filename',
				 'sim_epi_info_filename',
				 'read_sam_filename',
				 'sorted_bam_filename',
				 'annotated_info_filename',
				 'sim_bed_filename',
				 'sample',
				 'unique'
				 ]


	df_rows = []
	for exp in config:

		# generate conditions - every combination of desired parameters
		rows = list(itertools.product(
			range(config[exp]['replicates']), 	#2
			config[exp]['hosts'].keys(),	  	#3
			config[exp]['viruses'].keys(),		#4
			config[exp]['int_num'],
			config[exp]['epi_num'],
			config[exp]['min_sep'],
			config[exp]['min_len'],
			config[exp]['p_whole'],
			config[exp]['p_rearrange'],
			config[exp]['p_delete'],
			config[exp]['lambda_split'],
			config[exp]['p_overlap'],
			config[exp]['p_gap'],
			config[exp]['lambda_junction'],
			config[exp]['p_host_deletion'],
			config[exp]['lambda_host_deletion'],
			config[exp]['read_len'],
			config[exp]['fcov'],
			config[exp]['frag_len'],
			config[exp]['frag_std'],
			config[exp]['seq_sys']
			))
	
		# get number of conditions by removing the first element from each list
		n_conditions = len(set([i[1:] for i in rows]))
		condition_names = [f"cond{n}" for n in range(n_conditions)]*config[exp]['replicates']

		# add experiment and condition to conditions
		rows = [[exp, cond] + list(row) for cond, row in zip(condition_names, rows)]

		# add other derived information
		for i, row in enumerate(rows):

			# outpath
			row.append(config[exp]['out_directory'])
	
			# host_fasta
			row.append(ref_dict[row[3]])
	
			# virus_fasta
			row.append(ref_dict[row[4]])
	
			# random_seed
			row.append(config[exp]['initial_seed'] + i * config[exp]['seed_increment'])
			
			# genome simulation parameters
			row.append(config[exp]['snp_count'])		
			row.append(config[exp]['indel_count'])	
			row.append(config[exp]['cnv_count'])	
			row.append(config[exp]['inversion_count'])	
			row.append(config[exp]['translocation_count'])
	
			# sim_fa_filename
			row.append(f"{config[exp]['out_directory']}/{exp}/sim_ints/{row[1]}.rep{row[2]}.fa")
	
			# sim_int_info_filename
			row.append(f"{config[exp]['out_directory']}/{exp}/sim_ints/{row[1]}.rep{row[2]}.int-info.tsv")
	
			# sim_epi_info_filename
			row.append(f"{config[exp]['out_directory']}/{exp}/sim_ints/{row[1]}.rep{row[2]}.epi-info.tsv")
	
			# read_sam_filename
			row.append(f"{config[exp]['out_directory']}/{exp}/sim_reads/{row[1]}.rep{row[2]}.sam")
			
			# sorted_bam_filename
			row.append(f"{config[exp]['out_directory']}/{exp}/sim_reads/{row[1]}.rep{row[2]}.sorted.bam")
	
			# annotated_info_filename
			row.append(f"{config[exp]['out_directory']}/{exp}/sim_ints/{row[1]}.rep{row[2]}.int-info.annotated.tsv")
			
			# sim bed filename
			row.append(f"{config[exp]['out_directory']}/{exp}/sim_ints/{row[1]}.rep{row[2]}.int-info.bed")			
			
			# sample name: combination of condition and replicate
			sample = f"{row[1]}.rep{row[2]}"
			row.append(sample)
			
			# unique identifier
			row.append(f"{row[0]}__{sample}")

	
			df_rows.append(row)


	# create pandas dataframe from all rows
	df = pd.DataFrame(df_rows, columns = col_names)	

	return config, df, ref_dict
	
#### checking functions ####

# check that required field exists
def check_required(config, exp, field):
	if field not in config[exp]:
		raise ValueError(f"Please specify {field} for experiment {exp} in config file")

# check type of a field
def check_type(config, exp, field, type):
	if not(isinstance(config[exp][field], type)):
		raise ValueError(f"Field {field} in experiment {exp} must be of type {type}: (it's currently {type(config[exp][field])})")

# check that a field has a length of at least 1	
def check_at_least_one(config, exp, field):
	if len(config[exp][field]) < 1:
		raise ValueError(f"field {field} in experiment {exp} should have at least one entry")


# check if a field is set, otherwise set it with a default
def check_and_set_default(config, exp, field, default):
	if field not in config[exp]:
		config[exp][field] = default
		print(f"{field} not set in experiment {exp}, using {default}")
	return config

# check each entry in the list is of specified type
def check_each_type(config, exp, field, type):
	if not(all([isinstance(i, type) for i in config[exp][field]])):
		raise ValueError(f"{field} in experiment {exp} should have entries only of type {type}")

# make a list contain only unique values, if it does not already
def make_unique(config, exp, field):
	# if we have some redundant values
	if len(config[exp][field]) != len(set(config[exp][field])):
		config[exp][field] = list(set(config[exp][set]))
	return config

# check that an entry exists, has a specified type, and has at least one entry
def check_exists_type_length(config, exp, field, type):
	check_required(config, exp, field)
	check_type(config, exp, field, type)
	check_at_least_one(config, exp, field)


# check that entry exists, has a specified type, has at least one entry, and all entries are of a specified type
# if the entry does not exist, set it to a default
def standard_checks(config, exp, field, default, field_type, entry_type):
	config = check_and_set_default(config, exp, field, default)
	check_required(config, exp, field)
	check_type(config, exp, field, field_type)
	check_at_least_one(config, exp, field)
	check_each_type(config, exp, field, entry_type)
	config = make_unique(config, exp, field)
	
	return config

def standard_checks_int(config, exp, field, default):
		config = check_and_set_default(config, exp, field, default_replicates)
		check_required(config, exp, field)
		check_type(config, exp, field, int)
		
		return config
