
rule write_summary:
	input:
		list(ref_dict.values())
	output:
		tsv = "{outpath}/{exp}/simulation_summary.tsv"
	run:
		df[df['experiment'] == wildcards.exp].to_csv(output.tsv, sep='\t', index=False)


def get_parameter(wildcards, prefix, column_name):
	unique = f"{wildcards.exp}__{wildcards.cond}__{wildcards.rep}"
	value = df.loc[(df['unique'] == unique).idxmax(), column_name]
	return f"{prefix} {value}"

rule simulate_integrations:
	input:  
		host = lambda wildcards: df.loc[(df['unique'] == f"{wildcards.exp}__{wildcards.cond}__{wildcards.rep}").idxmax(), 'host_fasta'],
		virus =  lambda wildcards: df.loc[(df['unique'] == f"{wildcards.exp}__{wildcards.cond}__{wildcards.rep}").idxmax(), 'virus_fasta']
	output:
		sim_fasta = "{outpath}/{exp}/sim_ints/{cond}.rep{rep}.fa",
		sim_info = "{outpath}/{exp}/sim_ints/{cond}.rep{rep}.int-info.tsv",
		epi_info = "{outpath}/{exp}/sim_ints/{cond}.rep{rep}.epi-info.tsv",
	conda:
		"../envs/simvi.yml"
	container:
		"docker://szsctt/simvi:2"
	resources:
		mem_mb= lambda wildcards, attempt: attempt * 5000
	params:
		int_num = lambda wildcards: get_parameter(wildcards, '--int_num', 'int_num'),
		epi_num = lambda wildcards: get_parameter(wildcards, '--epi_num', 'epi_num'),
		p_whole = lambda wildcards: get_parameter(wildcards, '--p_whole', 'p_whole'),
		p_rearrange = lambda wildcards: get_parameter(wildcards, '--p_rearrange', 'p_rearrange'),
		p_delete = lambda wildcards: get_parameter(wildcards, '--p_delete', 'p_delete'),
		lambda_split = lambda wildcards: get_parameter(wildcards, '--lambda_split', 'lambda_split'),
		p_overlap = lambda wildcards: get_parameter(wildcards, '--p_overlap', 'p_overlap'),
		p_gap = lambda wildcards: get_parameter(wildcards, '--p_gap', 'p_gap'),
		lambda_junction = lambda wildcards: get_parameter(wildcards, '--lambda_junction', 'lambda_junction'),
		p_host_deletion = lambda wildcards: get_parameter(wildcards, '--p_host_deletion', 'p_host_deletion'),
		lambda_host_deletion = lambda wildcards: get_parameter(wildcards, '--lambda_host_deletion', 'lambda_host_deletion'),
		seed = lambda wildcards: get_parameter(wildcards, '--seed', 'random_seed'),
		
	shell:
		"""
		python3 scripts/insert_virus.py \
		 --host {input.host} \
		 --virus {input.virus} \
		 --ints {output.sim_fasta} \
		 --int_info {output.sim_info} \
		 --epi_info {output.epi_info} \
		 {params} \
		 --verbose
		"""

