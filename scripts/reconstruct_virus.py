# compare expected sequence from a viral fasta and information about integrations
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
	parser.add_argument('--virus-fa', help='fasta file of virus before integration', required = True, type=str)
	args = parser.parse_args()
	
	# import virus and integrated fasta
	virus = SeqIO.index(args.virus_fa, 'fasta', alphabet=unambiguous_dna)
	ints = SeqIO.index(args.int_fa, 'fasta', alphabet=unambiguous_dna)
	
	# loop over rows in info file and compare bases
	with open(args.int_info, newline = '') as csvfile:
		reader = csv.DictReader(csvfile, delimiter = '\t')
		for row in reader:
		
			# get integrated viral bases
			int_start = int(row['leftStop'])
			int_stop = int(row['rightStart'])
			int_chr = row['chr']
			int_bases = str(ints[int_chr][int_start:int_stop].seq)
			
			# get a list of coordinates within the virus
			virus_name = row['virus']
			pieces = eval(row['vBreakpoints'])
			
			# if we only have one piece, type of pieces[0] is 'int'
			if isinstance(pieces[0], int):
			
				# get viral bases
				virus_start = pieces[0]
				virus_stop = pieces[1]
				
				if row['vOris'] == "+":
					virus_bases = str(virus[virus_name][virus_start:virus_stop].seq)
				elif row['vOris'] == "-":
					virus_bases = str(virus[virus_name][virus_start:virus_stop].seq.reverse_complement())
				else:
					raise ValueError(f"unrecgonised orientation in row {row['id']}: {row['vOris']}")
				
			# otherwise we have to get multiple pieces
			elif isinstance(pieces[0], list):
				
				# get bases from virus
				virus_bases = ""
				
				# get orientations
				oris = row['vOris'].split(",")
				
				for i, piece in enumerate(pieces):
					virus_start = piece[0]
					virus_stop = piece[1]
					
					if oris[i] == "+":
						virus_bases = virus_bases + str(virus[virus_name][virus_start:virus_stop].seq)
					elif oris[i] == "-":
						virus_bases = virus_bases + str(virus[virus_name][virus_start:virus_stop].seq.reverse_complement())
					else:
						raise ValueError(f"unrecgonised orientation in row {row['id']}: {oris[i]}")
					
			else:
				raise ValueError(f"Error extracting breakpoints from row {row}")
			
			# compare virus with integrated virus
			if virus_bases.lower() != int_bases.lower():
				print(f"error with integration id {row['id']}: virus bases from {virus_start} to {virus_stop} were {virus_bases}, integrated bases from {int_start} to {int_stop} were {int_bases}")
			
	print("finished checking viruses")		
			
	
if __name__ == "__main__":
	main(argv[1:])
	
