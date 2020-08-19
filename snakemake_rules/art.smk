from os import path


rule art:
	input:
		sim_fasta = rules.simulate_integrations.output.sim_fasta,
	output:
		sam = "{outpath}/{exp}/sim_reads/{samp}.sam"
	params:
		seq_sys = lambda wildcards: get_parameter(wildcards, '-ss', 'seq_sys'),
		read_len = lambda wildcards: get_parameter(wildcards, '-l', 'read_len'),
		fcov = lambda wildcards: get_parameter(wildcards, '-f', 'fcov'),
		frag_len = lambda wildcards: get_parameter(wildcards, '-m', 'frag_len'),
		frag_std = lambda wildcards: get_parameter(wildcards, '-s', 'frag_std'),
		seed = lambda wildcards: get_parameter(wildcards, '--rndSeed', 'random_seed'),
		paried = "-p",
		sam = "-sam",
		input = lambda wildcards, input: f"-i {input.sim_fasta}",
		output = lambda wildcards, output: f"-o {path.splitext(output.sam)[0]}"
	conda:
		"../envs/art.yml"
	container:
		"docker://szsctt/art:1"
	resources:
		mem_mb= lambda wildcards, attempt: attempt * 10000,
		time = "24:00:00",
		nodes = 1
	shell:
		"""
		art_illumina {params}
		"""
		
rule convert:
	input:
		sam = rules.art.output.sam
	output:
		bam = "{outpath}/{exp}/sim_reads/{samp}.sorted.bam",
		idx = "{outpath}/{exp}/sim_reads/{samp}.sorted.bam.bai"
	conda:
		"../envs/bwa.yml"
	container:
		"docker://szsctt/bwa:1"
	resources:
		mem_mb= lambda wildcards, attempt: attempt * 5000,
		time = "24:00:00",
		nodes = 1
	shell:
		"""
		samtools sort -o {output.bam} {input.sam}
		samtools index {output.bam}
		"""
