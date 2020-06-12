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
import csv

###
max_attempts = 50 #maximum number of times to try to place an integration site 
int_virus_portions = ['whole', 'portion'] # possible portions of the virus to be integrated
int_junc_types = ['clean', 'gap']
int_junc_random_types = ['clean', 'gap'] # if an integration is one of these types, the position in the chromosome is random


def main(argv):
	#get arguments
	parser = argparse.ArgumentParser(description='simulate viral integrations')
	parser.add_argument('--host', help='host fasta file', required = True, type=str)
	parser.add_argument('--virus', help = 'virus fasta file', required = True, type=str)
	parser.add_argument('--ints', help = 'output fasta file', type=str, default="integrations.fa")
	parser.add_argument('--info', help = 'output tsv with info about viral integrations', type=str, default="integrations.tsv")
	parser.add_argument('--int_num', help = 'number of integrations to be carried out [5]', type=int, default=5)
	#parser.add_argument('--p_whole', help = 'probability of a virus being integrated completely', type = float, default=0.3) #TODO
	#parser.add_argument('--p_rearrange', help='probability of an integrated piece of virus being rearranged', type=float, default=0.01) #TODO
	#parser.add_argument('--p_delete', help='probability of an integrated piece of virus containing a deletion', type=float, default=0.01) #TODO
	parser.add_argument('--seed', help = 'seed for random number generator', default=1, type=int)
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
	if args.verbose is True:
		print("initialising a new Integrations object")
	ints = Integrations(host, virus, int_virus_portions)
	
	# do desired number of integrations
	if args.verbose is True:
		print(f"performing {args.int_num} integrations")
	for i in range(args.int_num):
		ints.add_integration(type = "rand")
	
	print("")
	if args.verbose is True:
		print(f"saving host fasta with integrations to {args.ints}")
	ints.save_fasta(args.ints)
	
	if args.verbose is True:
		print(f"saving information about integrations to {args.info}")
	ints.save_info(args.info)
	

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
		integration = Integration(self.host, self.virus, type)
		
		# append to self if nothing went wrong with this integration
		if integration.chunk.pieces is not None:
			if integration.pos is not None:
				self.append(integration)
		

		
	def save_fasta(self, filename):
		"""
		Save host fasta with integrated viral bases to a fasta file
		"""
		with open(filename, "w+") as handle:
			# loop over chromosomes
			for chr in self.host.keys():
			
				# print chromosome name
				handle.write(f">{chr}\n")
				position = 0
				
				# get integrations on this chromosome
				ints_on_this_chr = [integration for integration in self if integration.chr == chr]
				
				# if there are no integrations on this chromosome, just write the host bases
				if len(ints_on_this_chr) == 0:
					handle.write(str(self.host[chr].seq))
				
				# loop over sorted integrations in this chromosome
				else:
					for integration in sorted(ints_on_this_chr):
						# write host bases before this integration
						handle.write(str(self.host[chr].seq[position:integration.pos]))
						# write integrated viral bases
						handle.write(integration.chunk.bases)
						# update positon
						position = integration.pos
				# write the final bases of the chromosome
				handle.write(str(self.host[chr].seq[position:]))
				handle.write("\n")
		
	def save_info(self, filename):
		"""
		Output the following info for each integration (one integration per line) into a tab-separated file:
		Note that all co-orindates are 0-based
		 - chr: host chromosome on which integration occurs
		 - hPos: position in the original host fasta at which integration occurs
		 - hStart: start position of viral bases in host chromosome, accounting for any previous integrations
		 - hStop: stop position of viral bases in the host chromosome, accounting for any previous integrations
		 - virus: name of integrated virus
		 - viral breakpoints: a comma separated list of viral breakpoints which together indicate the integrated bases
			adjacent pairs of breakpoints indicate portions of the virus that have been integrated
		 - vOris: orientation (+ or -) of each portion of the virus that was integrated
		"""
		# dictionary will keep track of the number of bases previously integrated on each chromosome
		previous_ints = {key:0 for key in self.host.keys()}
		
		self.sort()
		with open(filename, 'w', newline='') as csvfile:
			intwriter = csv.writer(csvfile, delimiter = '\t')
			intwriter.writerow(['chr', 'hPos', 'intStart', 'intStop', 'virus', 'vBreakpoints', 'vOris'])
			for i, integration in enumerate(self):
				
				# calculate start and stop position for this integration
				hStart = integration.pos + previous_ints[integration.chr] + 1
				hStop = hStart + len(integration.chunk.bases) - 1
				
				# update previous_ints
				previous_ints[integration.chr] += len(integration.chunk.bases)

				breakpoints = ",".join([str(i) for i in integration.chunk.pieces])
				oris = ",".join(integration.chunk.oris)
				intwriter.writerow([integration.chr, 
									integration.pos, 
									hStart, 
									hStop, 
									integration.chunk.virus,
									breakpoints,
									oris])
	
	def __str__(self):
		return f"Viral integrations object with {len(self)} integrations of viral sequences {list(self.virus.keys())} into host chromosomes {list(self.host.keys())}"

	def __repr__(self):
		return f"Object of type Integrations with properties {self}"
		

	
class Integration:
	"""
	Class to store the properties of an individual integration
	"""
	def __init__(self, host, virus, portionType = 'whole', overlapTypes = ('clean', 'clean')):
		"""
		Function to initialise Integration object
		portionType is 'whole' or 'portion' - the part of the virus that has been inserted
		overlapType is two-member tuple of 'clean', 'gap' or 'overlap' - defines the junction at each end of the integration
		"""
		# check inputs
		assert isinstance(virus, dict)
		assert portionType in int_virus_portions
	
		# get random chromosome on which to do insertion
		self.chr = np.random.choice(list(host.keys()))
		
		# get viral chunk
		assert portionType in int_virus_portions
		chunk = ViralChunk(virus, portionType)
		
		# chunk.virus is assigned None if something went wrong with the initialisation
		if chunk.virus is not None:
			self.chunk = ViralChunk(virus, portionType)
		else:
			self.pos = None
			return
		assert self.chunk.pieces is not None
		assert self.chunk.type == portionType
		
		# assign portionType
		self.portionType = portionType
		
		# set overlap types
		assert isinstance(overlapTypes, tuple)
		assert len(overlapTypes) == 2
		assert all([True if (i in int_junc_types) else False for i in overlapTypes])
		self.overlapTypes = overlapTypes
		
		# get a position at which to integrate
		self.get_int_position(len(host[self.chr].seq), self.overlapTypes[0], self.overlapTypes[1])
		
	def get_int_position(self, chr_len, left_overlap, right_overlap,):
		"""
		get a position at which to perform the integration
		"""
		# if at both ends the junction is either 'clean' or 'gap', just get a random positon
		if all([True if (i in int_junc_random_types) else False for i in self.overlapTypes]):
			self.pos = np.random.randint(chr_len)
		
		
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
		assert isinstance(self.pos, int) and isinstance(other.pos, int)
		
		return (self.chr.lower(), self.pos) < (other.chr.lower(), other.pos)
		
	def __eq__(self, other):
		"""
		Integrations can be ranked by both chromosome name and position on that chromosome
		"""
		# make sure that we're comparing with another integration
		assert isinstance(other, Integration)
		
		# first check chromosome name
		assert isinstance(self.chr, str) and isinstance(other.chr, str)
		assert isinstance(self.pos, int) and isinstance(other.pos, int)
		
		return (self.chr.lower(), self.pos) == (other.chr.lower(), other.pos)

	
class ViralChunk:
	"""
	Class to store properties of an integrated chunk of virus
	"""
	
	def __init__(self, virus, type, deletion=True, rearrange=True, lambda_poisson = 1.5):
		"""
		function to get a chunk of a virus
		virus is the dictionary of seqrecords imported using biopython
		type specifies if we want the 'whole' virus, or just a 'portion'
		
		a viral chunk can be subdivided into smaller pieces if we want
		a deletion or a rearrangement.  if deletion is true, the chunk
		will be subdivided into at least three pieces and one of the middle ones
		will be deleted
		if rearrange is True, the chunk will be subdivided into at least two pieces
		and they will be randomly shuffled

		
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
			
		# breakpoints are the start and stop coordinates of pieces of the virus that have been 
		# integrated
		
		# set breakpoints
		self.pieces = ((self.start, self.stop),)
		
		self.oris = np.random.choice(('+', '-'))
		
		# do a deletion if applicable
		if deletion is True:
			self.__delete(lambda_poisson)
			# if something went wrong, breakpoints will be None
			if self.pieces is None:
				return
				
		if rearrange is True:
			self.__rearrange(lambda_poisson)
			# if something went wrong, breakpoints will be None
			if self.pieces is None:
				return
		
		# get bases
		self.bases = self.__get_bases(virus)
				
	def __get_bases(self, virus):
		"""
		return bases of viral chunk as a string
		note that objects of class ViralChunk do not store the whole viral sequence, so 
		need to provide this in order to get the bases
		"""
		bases = []
		# get integrated bases
		for i, points in enumerate(self.pieces):
			start = points[0]
			stop = points[1]
			if self.oris[i] == '+':
				bases += [str(virus[self.virus].seq[start:stop])]
			else:
				bases += [str(virus[self.virus].seq[start:stop].reverse_complement())]
		return "".join(bases)
		
	def __get_poisson_with_minimum_and_maximum(self, lambda_poisson, minimum = -np.inf, maximum = np.inf):
		"""
		A recurring task is to get an integer from the poisson distriubtion that is at least
		a minimum number
		"""
		# make sure minimum is less than or equal to maximum
		assert maximum >= minimum
		
		# if minimum and maximum are the same, just return 
		if minimum == maximum:
			return minimum
		
		# get number of pieces to divide chunk into
		counter = 0
		assert lambda_poisson > 1
		while True:
			# get a number from the poisson distrubtion
			n = np.random.poisson(lambda_poisson)
			# if it's more than the minimum, we're done
			if n >= minimum and n <= maximum:
				return n
			counter += 1
			# after too many tries, give up
			if counter > max_attempts:
				self.pieces = None
				return
		# if something went wrong and we ended up here, just set self.pieces to None
		# and return
		self.pieces = None
		return 0
		
	def __get_unique_list(self, start, stop, num = 1):
		"""
		get a random list of num integers between start and stop, with no repeats
		"""
		# make sure that we're not trying to get more numbers than is possible
		assert isinstance(start, int) and isinstance(stop, int) and isinstance(num, int)
		assert num <= (stop - start)
		
		return np.random.choice(range(start, stop), num, replace = False)
					
		
	def __split_into_pieces(self, lambda_poisson, min_pieces = 2):
		"""
		get random, unique breakpoints to divide a viral chunk into pieces
		there must be at least min_breakpoints, which results in min_breakpoints + 1 pieces 
		this is a list of coordinates, not tuples (unlike self.pieces)
		"""
		# shouldn't do this if this chunk has already been split
		assert len(self.pieces) == 1
		assert len(self.oris) == 1
		
		# check that we're not trying to divide a chunk into more pieces than there are bases
		if min_pieces >= self.stop - self.start:
			self.pieces = None
			return
		
		# get the number of pieces to divide into
		num_breakpoints = self.__get_poisson_with_minimum_and_maximum(lambda_poisson, minimum = min_pieces - 1)
		if num_breakpoints == 0:
			self.pieces = None
			return
				
		# get random breakpoints from within this chunk
		breakpoints = self.__get_unique_list(self.start + 1, self.stop - 1, num_breakpoints)
				
		# set self.pieces
		breakpoints = [self.start] + sorted(breakpoints) + [self.stop]
		self.pieces = [(breakpoints[i], breakpoints[i+1]) for i in range(len(breakpoints) - 1)]
		
		# set self.oris
		
		self.oris = [self.oris[0]] * len(self.pieces)
		return	
		
	def __swap_orientations(self, breakpoint, side = 'left'):
		"""
		Given a breakpoint, swap all of the orientations (+ to - or vice versa) for all of the pieces
		on the left or right of this breakpoint
		"""
		if side == 'left':
			for i, ori in self.oris[:breakpoint]:
				if ori == "+":
					self.oris[i] = "-"
				else:
					self.oris[i] = "+"
		else:
			for i, ori in self.oris[breakpoint:]:
				if ori == "+":
					self.oris[i] = "-"
				else:
					self.oris[i] = "+"
			
			
	
	def __delete(self, lambda_poisson):
		"""
		Divide a viral chunk up into multiple pieces
		and remove one of those pieces
		"""
		# deletions are always performed first, so the chunk should not have been split yet
		assert len(self.pieces) == 1

		# split chunk into at least three pieces
		self.__split_into_pieces(lambda_poisson, min_pieces = 3)
		if self.pieces is None:
			return
			
		# decide how many portions to delete - the maximum is the length of self.pieces - 2
		max_delete = len(self.pieces) - 2
		n_delete = self.__get_poisson_with_minimum_and_maximum(lambda_poisson, minimum = 0, maximum = max_delete)
		
		# decide which portions to delete
		i_delete = self.__get_unique_list(1, len(self.pieces) - 1, n_delete)
		
		# do deletion
		self.pieces = [piece for i, piece in enumerate(self.pieces) if i not in i_delete]
		self.oris = [ori for i, ori in enumerate(self.oris) if i not in i_delete]	
		
	def __rearrange(self, lambda_poisson):
		"""
		Divide a viral chunk up into multiple pieces
		and randomise their order and orientiations
		"""
		pdb.set_trace()
		# split the chunk if it hasn't already been split
		if len(self.pieces) == 1:
			# split chunk into at least three pieces
			self.__split_into_pieces(lambda_poisson, min_pieces = 2)
			if self.pieces is None:
				return
			
		# decide how many swaps to make
		n_swap = self.__get_poisson_with_minimum_and_maximum(lambda_poisson, minimum = 1)
		
		for i in range(n_swap):
			# pick a point about which to swap
			i_swap = self.__get_poisson_with_minimum_and_maximum(lambda_poisson, minimum = 1, maximum = len(self.pieces) - 1)
			
			# swap everything to the left of this position with everything on the right
			self.pieces = self.pieces[i_swap:] + self.pieces[:i_swap]
			
			# 50 % chance of swapping the orientations of all the pieces for each side
			if np.random.choice((True, False)) is True:
				self.__swap_orientations(self, i_swap, side = 'left')
			if np.random.choice((True, False)) is True:
				self.__swap_orientations(self, i_swap, side = 'right')

		
	
	def __str__(self):
		return f"Viral chunk of virus {self.virus} ({self.start}, {self.stop}) with type '{self.type}' and orientations {self.oris}"
		
	def __repr__(self):
		return f"Object of type ViralChunk with properties {self}"

if __name__ == "__main__":
	main(sys.argv[1:])
