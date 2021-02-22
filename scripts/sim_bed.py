#!/usr/bin/env python3

## make bed file with location of simulated virus/host junctions
## use output locations from simulate_integrations.py

## usage: python3 --sim_info <sim_info_file> --output <output_bed_file>

import sys
import argparse
import csv

def main(args):
	parser = argparse.ArgumentParser(description='Convert tab-separated file with integration locations into bed file with junction locations')
	parser.add_argument('--sim_info', '-i', help='Tab-separated output file with location of integrations', required=True)
	parser.add_argument('--output', '-o', help='Output bed file', default='ints.bed')
	args = parser.parse_args(args[1:])
	
	with open(args.sim_info, 'r', newline='') as infile, open(args.output, 'w', newline='') as outfile:
		inreader = csv.DictReader(infile, delimiter='\t')
		outwriter = csv.writer(outfile, delimiter='\t')
		
		for line in inreader:
		
			chrom = line['chr']
			pos = int(line['hPos'])

				
			# left ambiguous bases
			left_type = line['juncTypes'].split(',')[0]
			left_len = int(line['juncLengths'].split(',')[0])
			
			# check if left junction has supporting reads
			if not (line['left_chimeric'] == '' and line['left_discord'] == ''):			
				# write left junction
				#outwriter.writerow((chrom, pos, pos+left_len, '+'))	
				outwriter.writerow((chrom, pos, pos+left_len))

			# right ambiguous bases
			right_type = line['juncTypes'].split(',')[1]
			right_len = int(line['juncLengths'].split(',')[1])
			
			# deleted bases
			deleted = int(line['hDeleted'])
			
			# check if right junction has supporting reads			
			if not (line['right_chimeric'] == '' and line['right_discord'] == ''):			
				# write right junction
				#outwriter.writerow((chrom, pos+left_len+deleted, pos+left_len+deleted+right_len, '-'))
				outwriter.writerow((chrom, pos+left_len+deleted, pos+left_len+deleted+right_len))
	
if __name__ == "__main__":
	main(sys.argv)
