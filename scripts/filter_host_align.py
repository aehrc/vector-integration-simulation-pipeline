#!/usr/bin/env python3

# Filter host alignment (with all reads) for only those reads indicated to be integrations


from sys import argv
from os import path
import argparse
import csv
import pysam
import pdb

def main(argv):
	#get arguments
	parser = argparse.ArgumentParser(description='simulate viral integrations')
	parser.add_argument('--analysis-info', help='information from analysis of simulated reads', required = True, type=str)
	parser.add_argument('--sim-bam', help='alignment file with reads aligned to host', required = True, type=str)
	parser.add_argument('--output-bam', help='filtered bam', required=False)
	parser.add_argument('--chimeric', help='include chimeric reads', action='store_true')
	parser.add_argument('--discordant', help='include discordant read-pairs', action='store_true')
	args = parser.parse_args() 

	# default name for output file is input file + 'filtered'
	if args.output_bam is None:
		args.output_bam = path.os.splitext(args.sim_bam)[0] + 'filtered.bam'
		
	# make a set reads which are integrations
	ints = set()
	with open(args.analysis_info, newline="") as analysis:
		reader = csv.DictReader(analysis, delimiter='\t')
		for row in reader:
			ints.add(row['ReadID'])
			
	# filter alignment
	in_file = pysam.AlignmentFile(args.sim_bam)
	out_file = pysam.AlignmentFile(args.output_bam, 'wb', template=in_file)
	for read in in_file:
		if args.discordant is True:
			if read.qname in ints:
				out_file.write(read)
		if args.chimeric is True:
			if read.is_read1:
				if f"{read.qname}/1" in ints:
					out_file.write(read)
			elif read.is_read2:
				if f"{read.qname}/2" in ints:
					out_file.write(read)
			
			
	in_file.close()
	out_file.close()
	
	
	
if __name__ == "__main__":
	main(argv[1:])
