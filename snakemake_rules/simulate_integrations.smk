
rule write_summary:
	input:
		list(ref_dict.values())
	output:
		tsv = "{outpath}/{exp}/simulation_summary.tsv"
	run:
		df[df['experiment'] == wildcards.exp].to_csv(output.tsv, sep='\t', index=False)

def get_parameter(wildcards, column_name):
	unique = f"{wildcards.exp}__{wildcards.samp}"
	return df.loc[(df['unique'] == unique).idxmax(), column_name]


def format_parameter(wildcards, prefix, column_name):
	return f"{prefix} {get_parameter(wildcards, column_name)}"

rule simulate_integrations:
	input:  
		host = lambda wildcards: get_parameter(wildcards, 'host_fasta'),
		virus =  lambda wildcards:  get_parameter(wildcards, 'virus_fasta')
	output:
		sim_fasta = "{outpath}/{exp}/sim_ints/{samp}.fa",
		sim_info = "{outpath}/{exp}/sim_ints/{samp}.int-info.tsv",
		epi_info = "{outpath}/{exp}/sim_ints/{samp}.epi-info.tsv",
	conda:
		"../envs/simvi.yml"
	container:
		"docker://szsctt/simvi:2"
	resources:
		mem_mb= lambda wildcards, attempt, input: attempt * ( 5 * int((path.getsize(input.host)/1000000)) + int(get_parameter(wildcards, 'int_num')) * 50 + int(get_parameter(wildcards, 'epi_num')) * 10),
		time = lambda wildcards, attempt: ('2:00:00', '24:00:00', '24:00:00', '7-00:00:00')[attempt - 1],
		nodes = 1
	params:
		int_num = lambda wildcards: format_parameter(wildcards, '--int_num', 'int_num'),
		epi_num = lambda wildcards: format_parameter(wildcards, '--epi_num', 'epi_num'),
		p_whole = lambda wildcards: format_parameter(wildcards, '--p_whole', 'p_whole'),
		p_rearrange = lambda wildcards: format_parameter(wildcards, '--p_rearrange', 'p_rearrange'),
		p_delete = lambda wildcards: format_parameter(wildcards, '--p_delete', 'p_delete'),
		lambda_split = lambda wildcards: format_parameter(wildcards, '--lambda_split', 'lambda_split'),
		p_overlap = lambda wildcards: format_parameter(wildcards, '--p_overlap', 'p_overlap'),
		p_gap = lambda wildcards: format_parameter(wildcards, '--p_gap', 'p_gap'),
		lambda_junction = lambda wildcards: format_parameter(wildcards, '--lambda_junction', 'lambda_junction'),
		p_host_deletion = lambda wildcards: format_parameter(wildcards, '--p_host_deletion', 'p_host_deletion'),
		lambda_host_deletion = lambda wildcards: format_parameter(wildcards, '--lambda_host_deletion', 'lambda_host_deletion'),
		seed = lambda wildcards: format_parameter(wildcards, '--seed', 'random_seed'),
		
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

