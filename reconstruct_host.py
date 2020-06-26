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
			#print(f"checking row {row}")
			
			## if we're on a new chromosome
			if row['chr'] != current_chr:
				host_pos = 0
				int_pos = 0
				prev_right_start = 0
				prev_right_stop = 0
				prev_right_overlap = 0
				prev_deleted = 0
				current_chr = row['chr']
				
			## get left and right overlap types
			left_overlap = row['juncTypes'].split(",")[0] == "overlap"
			right_overlap = row['juncTypes'].split(",")[1] == "overlap"
		
			## get start and stop positions within integrated fasta
			int_start = int_pos
			
			# if we have a left overlap, we need to include the junction bases
			if left_overlap is True:
				int_stop = int(row['leftStop'])
			else:
				int_stop = int(row['leftStart'])
		
		
			## get start and stop positions within host
			
			host_start = host_pos + prev_deleted
			if left_overlap is True and prev_right_overlap is True:
				host_stop = host_start + (int(row['leftStop']) - prev_right_start)
			elif left_overlap is True:
				host_stop = host_start + (int(row['leftStop']) - prev_right_stop)
			elif prev_right_overlap is True:
				host_stop = host_start + (int(row['leftStart']) - prev_right_start)
			else:
				host_stop = host_start + (int(row['leftStart']) - prev_right_stop)
				
			
			# compare bases
			host_bases =  str(host[current_chr][host_start:host_stop].seq).lower()
			int_bases = str(ints[current_chr][int_start:int_stop].seq).lower()
			if host_bases != int_bases:
				print(f"error with integration {row['id']}: host bases from {host_start} to {host_stop} were {host_bases}, integrated bases from {int_start} to {int_stop} were {int_bases}")
			
			
			# update current positions
			if right_overlap is True:
				int_pos = int(row['rightStart'])
			else:
				int_pos = int(row['rightStop'])
				
			host_pos = host_stop - host_start
				
			prev_right_start = int(row['rightStart'])
			prev_right_stop = int(row['rightStop'])
			prev_right_overlap = right_overlap
			prev_deleted = int(row['hDeleted'])
			
			
	
if __name__ == "__main__":
	main(argv[1:])
	
