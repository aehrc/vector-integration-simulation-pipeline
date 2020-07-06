# compare expected sequence from a host fasta and information about integrations
# with the actual output


import csv
from Bio import SeqIO
from Bio.Alphabet.IUPAC import unambiguous_dna
from Bio.Seq import Seq 
from Bio.SeqRecord import SeqRecord
from sys import argv
import argparse
import pdb
import re

def main(argv):
	#get arguments
	parser = argparse.ArgumentParser(description='reconstruct host fasta after ')
	parser.add_argument('--int-fa', help='fasta file of host with integrations', required = True, type=str)
	parser.add_argument('--int-info', help = 'information file from integration script', required = True, type=str)
	parser.add_argument('--host-fa', help='fasta file of host before integration', required = True, type=str)
	args = parser.parse_args()
	
	# import host and integrated fasta
	host = SeqIO.index(args.host_fa, 'fasta', alphabet=unambiguous_dna)
	ints = SeqIO.index(args.int_fa, 'fasta', alphabet=unambiguous_dna)
	
	# loop over integrations and compare integrated fasta with host fasta
	current_chr = list(host.keys())[0]
	host_pos = 0
	int_pos = 0
	prev_right_start = 0
	prev_right_stop = 0
	prev_right_overlap = False
	prev_deleted = 0
	
	# loop over rows in info file and compare bases
	with open(args.int_info, newline = '') as csvfile:
		reader = csv.DictReader(csvfile, delimiter = '\t')
		for row in reader:
			
			## if we're on a new chromosome
			if row['chr'] != current_chr:
			
				# check between then end of the last integration and the end of the chromosome
				int_start = prev_right_stop
				int_stop = len(ints[current_chr])
				
				host_start = host_pos
				host_stop = len(host[current_chr])
				compare_bases(host[current_chr], ints[current_chr], host_start, host_stop, int_start, int_stop, f"bases at the end of chromosome {current_chr} (after integration {row['id']})")
			
				# update for a new chromosome
				host_pos = 0
				current_chr = row['chr']
				prev_right_stop = 0
				
			
			## compare bases between the previous integration and this one
			int_start = prev_right_stop
			int_stop = int(row['leftStart'])
			
			host_start = host_pos
			host_stop = host_start + (int(row['leftStart']) - prev_right_stop)
			host_pos = host_stop
			compare_bases(host[current_chr], ints[current_chr], host_start, host_stop, int_start, int_stop,  f"bases to the left of integration with id {row['id']}")
			
			## compare left overlap bases, if relevant
			if row['juncTypes'].split(",")[0] == "overlap":	
				# get coordinates in integrated fasta
				int_start = int(row['leftStart'])
				int_stop = int(row['leftStop'])
			
				# get coordinates in host fasta
				host_start = host_pos
				host_stop = host_start + int(row['juncLengths'].split(",")[0])
				host_pos = host_stop
			
				# compare integrated and host bases				
				compare_bases(host[current_chr], ints[current_chr], host_start, host_stop, int_start, int_stop, f"left overlap for integration with id {row['id']}")
		
			
			## compare right overlap bases, if relevant
			if row['juncTypes'].split(",")[1] == "overlap":
				# get coordinates in integrated fasta
				int_start = int(row['rightStart'])
				int_stop = int(row['rightStop'])
				
				# get coordinates in host fasta
				host_start = host_pos
				host_stop = host_start + int(row['juncLengths'].split(",")[1])
				host_pos = host_stop
			
				# compare integrated and host bases				
				compare_bases(host[current_chr], ints[current_chr], host_start, host_stop, int_start, int_stop, f"right overlap for integration with id {row['id']}")
		
			## keep track of the end of the this integration, for evaluating the next one
			host_pos = host_pos + int(row['hDeleted'])
			prev_right_stop = int(row['rightStop'])
			
	## compare the bases at the end of the last chromosome
	int_start = prev_right_stop
	int_stop = int(row['leftStart'])
			
	host_start = host_pos
	host_stop = host_start + (int(row['leftStart']) - prev_right_stop)
	compare_bases(host[current_chr], ints[current_chr], host_start, host_stop, int_start, int_stop, f"bases at the end of chromosome {current_chr} (after integration {row['id']})")

	print("finished checking host")		
			

def compare_bases(host_chr, int_chr, host_start, host_stop, int_start, int_stop, type = ''):
	"""
	compare host bases from host_start to host_stop to integrated bases from int_start to int_stop
	"""

	# get bases for host and integrated fasta chromosomes
	host_bases = str(host_chr[host_start:host_stop].seq).lower()
	int_bases = str(int_chr[int_start:int_stop].seq).lower()
	
	# compare
	if host_bases != int_bases:
		print(f"error with {type}: host bases from {host_start} to {host_stop} were {host_bases}, integrated bases from {int_start} to {int_stop} were {int_bases}")
		
	
if __name__ == "__main__":
	main(argv[1:])
	
