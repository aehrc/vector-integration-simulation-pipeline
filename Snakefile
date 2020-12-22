import sys
import os
import pdb

# parse config file and create data frame containing parameters for running analysis (from config file)

from snakemake_rules import parse_config

# apply global options
if 'global' in config:		
	# get default (global) options
	default = config.pop('global')
	for dataset in config:
		for key in default:
			if key not in config[dataset]:
				config[dataset][key] = default[key]

config, sim_df, ref_dict = parse_config(config)

#####################################################
################### target files ####################
#####################################################

rule all:
	input: 
		expand("{outpath}/{exp}/simulation_summary.tsv", 
			zip,
			exp = list(config.keys()), # config is imported as ordered dict
			outpath = [config[exp]['out_directory'] for exp in config.keys()]
			),
		sim_df.loc[:,'annotated_info_filename'],
		


#####################################################
############### simulate integrations ###############
#####################################################


include: "snakemake_rules/simulate_integrations.smk"

include: "snakemake_rules/art.smk"

include: "snakemake_rules/annotate_reads.smk"
