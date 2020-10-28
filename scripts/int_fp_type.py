#!/usr/bin/env python3

# integrations that are scored as false positives may consist of reads that are either:
# 1. do cross a host/virus junction, but are aligned in the wrong place, or
# 2. don't cross a host/virus junction
# try to distingish between the two by checking, for each of the reads in the integration
# whether or not the read was indicated to cross an integration

from sys import argv
import argparse
import csv
import re
import pdb


def main(argv):
	#get arguments
	parser = argparse.ArgumentParser(description='simulate viral integrations')
	parser.add_argument('--sim-info', help='information from simulation annotated with reads crossing junctions', required = True, type=str)
	parser.add_argument('--scored-ints', help='information from found of simulated reads', required = True, type=str)
	parser.add_argument('--merged', help='a list of merged readIDs, one per line', required = False, type=str)
	parser.add_argument('--output', help='output file', required=False, default='results.tsv')
	
	args = parser.parse_args()
	
	# for each row in sim, create a dict that contains all the reads (of any type) that cross either junction
	# of that integration
	sim = parse_sim(args.sim_info)
	
	# get all the merged readIDs, if provided
	merged_reads = set()
	if args.merged:
		with open(args.merged) as merged_file:
			for line in merged_file:
				merged_reads.add(line.strip())
	
	# read scored integrations and check false positives cross any junctions
	with open(args.scored_ints, newline = "") as infile, open(args.output, "w", newline = "") as outfile:
	
		# csv objects to read infile and write outfile
		reader = csv.DictReader(infile, delimiter = '\t')
		out_fieldnames = reader.fieldnames + ['n_junc_reads', 'n_non_junc_reads']
		writer = csv.DictWriter(outfile, delimiter = '\t', fieldnames = out_fieldnames)
		writer.writeheader()
		
		# find fp rows in scored ints
		for row in reader:

			if row['score'] != 'fp':
				row['n_junc_reads'] = ""
				row['n_non_junc_reads'] = ""
				writer.writerow(row)
				continue	
				
			# get reads for this integration
			reads = row['total_reads'].split(";")
			
			n_found = 0
			n_not_found = 0
			# check all the reads that are part of the fp integration
			for read in reads:
				found = False
				
				# was this read merged?
				if read in merged_reads:
					merged = True
				else:
					merged = False
					
			
				# for each simulated integration check if fp read crosses
				for sim_reads in sim.values():
				
					# if the read was merged, we need to strip '/1' and '/2' from reads
					if merged:
						# this removes any reads that don't have "/", but we don't care about these
						# since merged reads can only be chimeric, and only chmieric reads have "/"
						if read in {read[:-2] for read in sim_reads if read.find("/")}:
							found = True
							n_found += 1
							break
					else:	
						if read in sim_reads:
							n_found += 1
							found = True
							break
				
				# if we didn't find this read
				if not found:
					n_not_found += 1
					
			# output the number of found and not found reads
			row['n_junc_reads'] = n_found
			row['n_non_junc_reads'] = n_not_found
			writer.writerow(row)
			
	
	
	
	
def parse_sim(path):
	"""
	parse a tsv, and store each row as a dict in a list
	"""
	with open(path, newline = '') as handle:
		reader = csv.DictReader(handle, delimiter = '\t')
		data = {}
		
		# get reads and id from each row
		for row in reader:
			
			reads = ";".join((row['left_chimeric'], row['right_chimeric'], row['left_discord'], row['right_discord']))
			data[row['id']] = set(reads.split(';'))
			
	return data	

if __name__ == "__main__":
	main(argv)

