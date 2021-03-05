def return_lower(expr1, expr2):
	if expr1 <= expr2:
		return expr1
	return expr2

def bandpass(num, minimum, maximum):
	num = min(maximum, num)
	return max(num, minimum)

def get_mem(file_name_list, attempt, mult_factor=2, minimum=100, maximum=50000):
	resource = int(sum([os.stat(file).st_size/1e6 for file in file_name_list])) * mult_factor * attempt
	resource = min(maximum, resource)
	return max(minimum, resource)


rule write_summary:
	input:
		list(ref_dict.values())
	output:
		tsv = "{outpath}/{exp}/simulation_summary.tsv"
	run:
		sim_df[sim_df['experiment'] == wildcards.exp].to_csv(output.tsv, sep='\t', index=False)

def get_parameter(wildcards, column_name):
	unique = f"{wildcards.exp}__{wildcards.samp}"
	return sim_df.loc[(sim_df['unique'] == unique).idxmax(), column_name]

def format_parameter(wildcards, prefix, column_name):
	return f"{prefix} {get_parameter(wildcards, column_name)}"
	
def get_host_parameter(wildcards, column_name):
	exp_rows = sim_df.loc[(sim_df['experiment'] == wildcards.exp)]
	return exp_rows.loc[(exp_rows['host_name'] == wildcards.sim_host).idxmax(), column_name]

#rule genome_variation:
#	input:  
#		host = lambda wildcards: {name:fasta for name, fasta in zip(sim_df.host_name, sim_df.host_fasta)}[wildcards.sim_host]
#	output:
#		r1_sim_fasta = "{outpath}/references/{exp}/simuG/{sim_host}.1.simseq.genome.fa",
#		vcf_snp = temp("{outpath}/references/{exp}/simuG/{sim_host}.1.refseq2simseq.SNP.vcf"),
#		vcf_indel = temp("{outpath}/references/{exp}/simuG/{sim_host}.1.refseq2simseq.INDEL.vcf"),
#		r1_map = "{outpath}/references/{exp}/simuG/{sim_host}.1.refseq2simseq.map.txt"
#	container:
#		"docker://szsctt/simug:1"
#	resources:
#		mem_mb= lambda wildcards, attempt, input:get_mem(input, attempt, mult_factor=2, minimum=5000, maximum=50000),
#		time = lambda wildcards, attempt: ('30:00', '2:00:00', '24:00:00', '7-00:00:00')[attempt - 1],
#		nodes = 1
#	wildcard_constraints:
#		sim_host = "|".join(set(sim_df.loc[:, 'host_name'])),
#	params:
#		prefix = lambda wildcards, output: f"-prefix {os.path.splitext(os.path.splitext(os.path.splitext(output.r1_sim_fasta)[0])[0])[0]}",
#		snp_count = lambda wildcards: f"-snp_count {get_host_parameter(wildcards, 'snp_count')}",
#		indel_count = lambda wildcards: f"-indel_count {get_host_parameter(wildcards, 'indel_count')}",
#		#cnv_count = lambda wildcards: f"-cnv_count {get_host_parameter(wildcards, 'cnv_count')}",
#		#inversion_count = lambda wildcards: f"-inversion_count {get_host_parameter(wildcards, 'inversion_count')}",
#		#translocation_count = lambda wildcards: f"-translocation_count {get_host_parameter(wildcards, 'translocation_count')}",
#		seed = "-seed 202102231543"
#	shell:
#		"""
#		perl /var/work/simuG/simuG.pl -refseq {input.host} {params}
#		"""

rule simulate_integrations:
	input:  
		virus =  lambda wildcards: get_parameter(wildcards, 'virus_fasta'),
		host = lambda wildcards: get_parameter(wildcards, 'host_fasta')
#		host = lambda wildcards: expand(rules.genome_variation.output.r1_sim_fasta, sim_host = get_parameter(wildcards, 'host_name'), allow_missing=True),
#		r1_map = lambda wildcards: expand(rules.genome_variation.output.r1_map, sim_host = get_parameter(wildcards, 'host_name'), allow_missing=True),
	output:
		sim_fasta = temp("{outpath}/{exp}/sim_ints/{samp}.fa"),
		sim_info = temp("{outpath}/{exp}/sim_ints/{samp}.int-info.tsv"),
		epi_info = "{outpath}/{exp}/sim_ints/{samp}.epi-info.tsv",
	conda:
		"../envs/simvi.yml"
	container:
		"docker://szsctt/simvi:2"
	resources:
		mem_mb= lambda wildcards, attempt, input: bandpass(get_mem(input, attempt, 10) + int(get_parameter(wildcards, 'int_num'))*10 + int(get_parameter(wildcards, 'epi_num'))*10, 100, 50000),
		time = lambda wildcards, attempt: ('30:00', '2:00:00', '24:00:00', '7-00:00:00')[attempt - 1],
		nodes = 1
	params:
		int_num = lambda wildcards: format_parameter(wildcards, '--int_num', 'int_num'),
		epi_num = lambda wildcards: format_parameter(wildcards, '--epi_num', 'epi_num'),
		min_sep = lambda wildcards: format_parameter(wildcards, '--min-sep', 'min_sep'),
		min_len = lambda wildcards: format_parameter(wildcards, '--min-len', 'min_len'),
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
#		 --simug_snp_indel {input.r1_map} \
