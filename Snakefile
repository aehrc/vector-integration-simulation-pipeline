import sys
import os

# parse config file and create data frame containing parameters for running analysis (from config file)

from snakemake_rules import parse_config

config, df, ref_dict = parse_config(config)


#####################################################
################### target files ####################
#####################################################

rule all:
	input: 
		expand("{outpath}/{exp}/experiment_summary.tsv", 
			zip,
			exp = list(config.keys()), # config is imported as ordered dict
			outpath = [config[exp]['out_directory'] for exp in config.keys()]
			),
		df.loc[:,'read_sam_filename'],
		


#####################################################
############### simulate integrations ###############
#####################################################


include: "snakemake_rules/simulate_integrations.smk"
include: "snakemake_rules/art.smk"

