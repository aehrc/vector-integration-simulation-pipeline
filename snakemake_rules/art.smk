from os import path


rule art:
	input:
		sim_fasta = rules.simulate_integrations.output.sim_fasta,
	output:
		sam = temp("{outpath}/{exp}/sim_reads/{samp}.sam"),
		r1 = temp("{outpath}/{exp}/sim_reads/{samp}1.fq"),
		r2 = temp("{outpath}/{exp}/sim_reads/{samp}2.fq"),
		aln1 = temp("{outpath}/{exp}/sim_reads/{samp}1.aln"),
		aln2 = temp("{outpath}/{exp}/sim_reads/{samp}2.aln"),
	params:
		seq_sys = lambda wildcards: format_parameter(wildcards, '-ss', 'seq_sys'),
		read_len = lambda wildcards: format_parameter(wildcards, '-l', 'read_len'),
		fcov = lambda wildcards: format_parameter(wildcards, '-f', 'fcov'),
		frag_len = lambda wildcards: format_parameter(wildcards, '-m', 'frag_len'),
		frag_std = lambda wildcards: format_parameter(wildcards, '-s', 'frag_std'),
		seed = lambda wildcards: format_parameter(wildcards, '--rndSeed', 'random_seed'),
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
		time = lambda wildcards, attempt: ('30:00', '2:00:00', '24:00:00', '7-00:00:00')[attempt - 1],
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
		time = lambda wildcards, attempt: ('30:00', '2:00:00', '24:00:00', '7-00:00:00')[attempt - 1],
		nodes = 1
	shell:
		"""
		samtools sort -o {output.bam} {input.sam}
		samtools index {output.bam}
		"""
