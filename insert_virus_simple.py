###import libraries
from Bio import SeqIO
from Bio.Alphabet.IUPAC import unambiguous_dna
from Bio.Seq import Seq 
from Bio.SeqRecord import SeqRecord
import pandas as pd
import argparse
import sys
import os
import numpy as np
import pdb 

###
max_attempts = 50 #maximum number of times to try to place an integration site 
int_virus_portions = ['whole', 'portion'] # possible portions of the virus to be integrated
int_junc_types = ['clean']

def main(argv):
	#get arguments
	parser = argparse.ArgumentParser(description='simulate viral integrations')
	parser.add_argument('--host', help='host fasta file', required = True, type=str)
	parser.add_argument('--virus', help = 'virus fasta file', required = True, type=str)
	parser.add_argument('--ints', help = 'output fasta file', required = True, type=str)
	parser.add_argument('--locs', help = 'output text with viral integrations', required = True, type=str)
	parser.add_argument('--locs_host', help = 'output csv with integration locations in host genome', required = True, type=str)
	parser.add_argument('--int_num', help = 'number of integrations to be carried out [5]', required=True, type=int, default=5)
	parser.add_argument('--seed', help = 'seed for random number generator', required=False, default=1, type=int)
	parser.add_argument('--verbose', help = 'display extra output for debugging', action="store_true")
	args = parser.parse_args()
	
	
	if args.verbose is True:
		print("importing host fasta")
	
	#read host fasta - use index which doesn't load sequences into memory because host is large genome
	if checkFastaExists(args.host):
		host = SeqIO.index(args.host, 'fasta', alphabet=unambiguous_dna)
	else:
		raise OSError("Could not open host fasta")
		
	if args.verbose is True:
		print(f"imported host fasta with {len(host)} chromosomes:")
		host_chrs = {key:len(host[key]) for key in host.keys()}
		for key, length in host_chrs.items():
			print(f"\thost chromosome '{key}' with length {length}") 
	
	#read virus fasta -  make dictionary in memory of sequences
	if checkFastaExists(args.virus):
		virus = SeqIO.to_dict(SeqIO.parse(args.virus, 'fasta', alphabet=unambiguous_dna))
		# convert all viral sequences to upper case to make it easier to check output
		for this_virus in virus.values():
			this_virus.seq = this_virus.seq.upper()
	else:
		raise OSError("Could not open virus fasta")

	if args.verbose is True:
		print(f"imported virus fasta with {len(virus)} sequences:")
		virus_seqs = {key:len(virus[key]) for key in virus.keys()}
		for key, length in virus_seqs.items():
			print(f"\tviral sequence '{key}' with length {length}") 
			
	#set random seed
	np.random.seed(args.seed)
			
	# initialise integrations object
	ints = Integrations(host, virus, int_virus_portions)
	ints.add_integration()
	
	print("checking one integration")
	print(ints)
	print(ints[0])
	print(ints[0].chunk)
	print(len(ints))
	
	print("\nchecking two integrations")
	ints.add_integration()
	print(ints)
	print(ints[1])
	print(f"integration 0 more than integration 1: {ints[0] > ints[1]}")
	print(f"integration 1 more than integration 0: {ints[1] > ints[0]}")
	
	print("\nchecking sorting")
	ints.add_integration()
	ints.add_integration()
	ints.add_integration()
	ints.sort()
	[print(f"{ints[i]}") for i in range(len(ints))]
	
	print("\nchecking saving a fasta:")
	if args.verbose is True:
		print(f"saving host fasta with integrations to {args.ints}")
	ints.save_fasta(args.ints)
	
			

def checkFastaExists(file):
	#check file exists
	exists = os.path.isfile(file)
	if not(exists):
		return False
	#check extension
	prefix = file.split(".")[-1]
	if prefix:
		if not((prefix == "fa") or (prefix == "fasta") or (prefix == "fna")):
			return False
	return True
	
class Integrations(list):
	"""
	Class to store all integrations for a given host and virus
	"""
	def __init__(self, host, virus, int_types):
		"""
		Function to initialise Integrations object
		"""
		# assign properties common to all integrations
		self.host = host
		self.virus = virus
		self.int_types = int_types
		
	def add_integration(self, type = 'rand'):
		"""
		Add an integration
		"""
		# either get random type or check this this type is valid
		if type == 'rand':
			type = np.random.choice(self.int_types)
		else:
			assert type in self.types
			
		# make an integration
		int = Integration(self.host, self.virus, type)
		
		# append integration to self.ints
		self.append(int)
		
	def save_fasta(self, filename):
		"""
		Save all integrations to a fasta file
		"""
		
		with open(filename, "w+") as handle:
			# loop over chromosomes
			for chr in self.host.keys():
			
				# print chromosome name
				handle.write(f">{chr}\n")
				position = 0
				
				# loop over sorted integrations in this chromosome
				for int in sorted([int for int in self if int.chr == chr]):
					# write host bases before this integration
					handle.write(str(self.host[chr].seq[position:int.pos]))
					# write integrated viral bases
					handle.write(int.chunk.bases)
					# update positon
					position = int.pos
				handle.write("\n")
		
		
		
		
	def __str__(self):
		return f"Viral integrations object with {len(self)} integrations of viral sequences {list(self.virus.keys())} into host chromosomes {list(self.host.keys())}"

	def __repr__(self):
		return f"Object of type Integrations with properties {self}"
		

	
class Integration:
	"""
	Class to store the properties of an individual integration
	"""
	def __init__(self, host, virus, type):
		"""
		Function to initialise Integration object
		"""
		# check inputs
		assert isinstance(virus, dict)
		assert type in int_virus_portions
	
		# get random chromosome on which to do insertion
		self.chr = np.random.choice(list(host.keys()))
	
		# get random position at which to insert
		self.pos = np.random.randint(1, len(host[self.chr].seq))
		
		# assign type
		self.type = type
		
		# get viral chunk
		self.chunk = ViralChunk(virus, type)
		
		# check viral type
		assert self.chunk.type == self.type
		
		
	def __str__(self):
		return f"Viral integration into host chromosome {self.chr} at position {self.pos} of type '{self.type}'"
		
	def __repr__(self):
		return f"Object of type Integration with properties {self}"
		
	def __lt__(self, other):
	
		"""
		Integrations can be ranked by both chromosome name and position on that chromosome
		"""
		# make sure that we're comparing with another integration
		assert isinstance(other, Integration)
		
		# first check chromosome name
		assert isinstance(self.chr, str) and isinstance(other.chr, str)
		if self.chr > other.chr:
			return False
		
		# then check position
		assert isinstance(self.pos, int) and isinstance(other.pos, int)
		if self.pos > other.pos:
			return False
		
		# if not false, then true
		return True
	
class ViralChunk:
	"""
	Class to store properties of an integrated chunk of virus
	"""	
	
	def __init__(self, virus, type):
		"""
		function to get a chunk of a virus
		"""
		# check inputs
		assert isinstance(virus, dict)
		assert type in int_virus_portions
		
		# assign type to self
		self.type = type
		
		# get a random virus to integrate
		self.virus = np.random.choice(list(virus.keys()))
		
		if type == "whole":
			self.start = 0
			self.stop = len(virus[self.virus])
		else:
			self.start = np.random.randint(0, len(virus[self.virus].seq) -1)
			self.stop = np.random.randint(self.start + 1, len(virus[self.virus].seq))
			
		# get integrated bases
		self.bases = str(virus[self.virus].seq[self.start:self.stop])
			
			
	def __str__(self):
		return f"Viral chunk of virus {self.virus} ({self.start}, {self.stop}) with type '{self.type}'"
		
	def __repr__(self):
		return f"Object of type ViralChunk with properties {self}"

if __name__ == "__main__":
	main(sys.argv[1:])
