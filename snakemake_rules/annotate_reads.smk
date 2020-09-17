rule annotate_reads:
	input:
		bam = rules.convert.output.bam,
		info = rules.simulate_integrations.output.sim_info
	output:
		annotated_info = "{outpath}/{exp}/sim_ints/{samp}.int-info.annotated.tsv",
	conda:
		"../envs/simvi.yml"
	container:
		"docker://szsctt/simvi:2"
	resources:
		mem_mb= lambda wildcards, attempt: attempt * 10000,
		time = "24:00:00",
		nodes = 1
	params:
		mean_frag_len = lambda wildcards: format_parameter(wildcards, '--mean-frag-len', 'frag_len'),
		sd_frag_len = lambda wildcards: format_parameter(wildcards, '--sd-frag-len', 'frag_std'),
		window_frac = "--window-frac 0.99"
		
	shell:
		"""
		python3 scripts/annotate_reads.py \
		 --sim-info {input.info} \
		 --sim-bam {input.bam} \
		 --output {output.annotated_info} \
		 {params} 
		"""
