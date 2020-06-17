###import libraries
from Bio import SeqIO
from Bio.Alphabet.IUPAC import unambiguous_dna
from Bio.Seq import Seq 
from Bio.SeqRecord import SeqRecord
import pandas as pd
import argparse
from sys import argv
from os import path
import numpy as np
import pdb 
import csv
import re

###
max_attempts = 50 #maximum number of times to try to place an integration site 

default_ints = 5

default_p_whole = 0.3
default_p_rearrange = 0.05
default_p_delete = 0.05
default_lambda_split = 1.5
default_p_overlap = 0.3
default_p_gap = 0.3
default_lambda_junction = 3

search_length_overlap = 10000 # number of bases to search either side of randomly generated position for making overlaps


def main(argv):
	#get arguments
	parser = argparse.ArgumentParser(description='simulate viral integrations')
	parser.add_argument('--host', help='host fasta file', required = True, type=str)
	parser.add_argument('--virus', help = 'virus fasta file', required = True, type=str)
	parser.add_argument('--ints', help = 'output fasta file', type=str, default="integrations.fa")
	parser.add_argument('--info', help = 'output tsv with info about viral integrations', type=str, default="integrations.tsv")
	parser.add_argument('--int_num', help = f"number of integrations to be carried out [{default_ints}]", type=int, default=default_ints)
	parser.add_argument('--p_whole', help = f"probability of a virus being integrated completely [{default_p_whole}]", type = float, default=default_p_whole) 
	parser.add_argument('--p_rearrange', help=f"probability of an integrated piece of virus being rearranged [{default_p_rearrange}]", type=float, default=default_p_rearrange) 
	parser.add_argument('--p_delete', help=f"probability of an integrated piece of virus containing a deletion [{default_p_delete}]", type=float, default=default_p_delete) 
	parser.add_argument('--lambda_split', help = f"mean of poisson distriubtion for number of pieces to split an integrated viral chunk into when rearranging or deleting [{default_lambda_split}]", type=float, default=default_lambda_split)
	parser.add_argument('--p_overlap', help=f"probability of a junction to be overlapped (common bases between virus and host) [{default_p_overlap}]", type=float, default=default_p_overlap) 
	parser.add_argument('--p_gap', help=f"probability of a junction to have a gap [{default_p_gap}]", type=float, default=default_p_gap) 
	parser.add_argument('--lambda_junction', help = f"mean of possion distribution of number of bases in a gap/overlap [{default_lambda_junction}]", type=float, default=default_lambda_junction)
	parser.add_argument('--seed', help = 'seed for random number generator', default=12345, type=int)
	parser.add_argument('--verbose', help = 'display extra output for debugging', action="store_true")
	args = parser.parse_args()
	
	
	#### check and import fasta files ####
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
	
	#### generate integrations ####		
	
	# generate dictionary of probabilities and lambda/means as input for Integrations class
	probs = {'p_whole' 			: args.p_whole,
			'p_rearrange' 		: args.p_rearrange,
			'p_delete'			: args.p_delete,
			'lambda_split'		: args.lambda_split,
			'p_overlap' 		: args.p_overlap,
			'p_gap' 			: args.p_gap,
			'lambda_junction' 	: args.lambda_junction}
			
	# initialise integrations object
	if args.verbose is True:
		print("initialising a new Integrations object")
	ints = Integrations(host, virus, probs, seed=args.seed)
	
	# do desired number of integrations
	if args.verbose is True:
		print(f"performing {args.int_num} integrations")
	counter = 0
	while len(ints) < args.int_num:
		if ints.add_integration() is False:
			counter += 1
		if counter > max_attempts:
			raise valueError('too many failed attempts to add integrations')
	
	
	#### save outputs ####
	print("")
	if args.verbose is True:
		print(f"saving host fasta with integrations to {args.ints}")
	ints.save_fasta(args.ints)
	
	if args.verbose is True:
		print(f"saving information about integrations to {args.info}")
	ints.save_info(args.info)
	

def checkFastaExists(file):
	#check file exists
	exists = path.isfile(file)
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
	This class stores the sequences for the host and virus (dictionaries of biopython seqRecord objects)
	And probabilities of each type of viral integration - whole/portion, rearrangement/deletion, gap/overlap, 
	means of poisson distributions
	
	These probabilities and means are stored in a dictionary - probs, which must contain the following values:
	
	p_whole - probability of an integration being of the whole virus [0.3].  probability of a portion is 1-p_whole
	
	p_rearrange - probability of a chunk being rearranged
	p_delete - probability of a chunk containing a deletion
	
	note that rearrangements and deletions are not mutually exclusive - an integration can have both a rearrangement and 
	a deletion.  if this is the case, deletions are always performed first.
	
	lambda_split - when splitting a chunk for a rearrangement or deletion, mean of the poission distribution for
	number of pieces in which to split chunk
	
	p_overlap - probability of a junction to contain an overlap (common bases between virus and host)
	p_gap - probability of a junction to contain a gap (random bases inserted between host and virus)
	(probability of neither a rearrangement or a deletion is 1 - (p_overlap + p_gap) )
	lambda_junction - mean of poisson distribution of number of bases involved in overlap or gap
	
	"""
	def __init__(self, host, virus, probs, seed = 12345, max_attempts=50):
		"""
		Function to initialise Integrations object
		"""
		# assign properties common to all integrations
		self.host = host
		self.virus = virus
		self.probs = probs
		self.rng = np.random.default_rng(seed)
		self.max_attempts = max_attempts
		
		# default values for probs
		self.set_probs_default('p_whole', default_p_whole)
		self.set_probs_default('p_rearrange', default_p_rearrange)
		self.set_probs_default('p_delete', default_p_delete)
		self.set_probs_default('lambda_split', default_lambda_split)
		self.set_probs_default('p_overlap', default_p_overlap)
		self.set_probs_default('p_gap', default_p_gap)
		self.set_probs_default('lambda_junction', default_lambda_junction)		
			
		# check that probabilities are between 0 and 1
		self.check_prob(self.probs['p_whole'], 'p_whole')
		self.check_prob(self.probs['p_rearrange'] + self.probs['p_delete'], 'the sum of p_rearrange and p_delete')	
		self.check_prob(self.probs['p_overlap'] + self.probs['p_gap'], 'the sum of p_overlap and p_gap')
			
		# check that lambda values are positive floats	
		self.check_float_or_int_positive(self.probs['lambda_split'], 'lambda_split')
		self.check_float_or_int_positive(self.probs['lambda_junction'], 'lambda_junction')
	
	def set_probs_default(self, key, default):
		"""
		check if a key is in the probs dictionary, and if not set it to a default
		"""	
		if key not in self.probs:
			self.probs[key] = default
	
	def check_float_or_int_positive(self, value, name):
		"""
		Check that a value, with a name, is a positive float or int
		"""
		if not isinstance(value, float) and not isinstance(value, int):
			raise ValueError(f"{name} must be a float (it's currently {type(value)})")
		if value <= 0:
			raise ValueError(f"{name} must be greater than zero (it's currently {value})")	
			
	def check_prob(self, value, name):
		"""
		Check that a value, with a name, is a valid probability (float between 0 and 1):
		"""
		if not isinstance(value, float):
			raise ValueError(f"{name} must be a float (it's currently {type(value)})")
		if value < 0 or value > 1:
			raise ValueError(f"{name} must be between 0 and 1")	
			
	def poisson_with_minimum_and_maximum(self, lambda_poisson, min=-np.inf, max=np.inf):
		"""
		get an integer from the poisson distribution with the specified lambda
		with a minimum value 
		"""
		assert max >= min
			
		if min == max:
			return min
			
		counter = 0
		while True:
			n = int(self.rng.poisson(lam = lambda_poisson))
			if n >= min and n <= max:
				return n
			counter += 1
			if counter > self.max_attempts:
				return None
					
	def add_integration(self):
		"""
		Add an integration by appending an Integration object to self.  
		Decide the properties of the integration using the 'probs' dict
		"""
		
		# get if integration should be whole or portion
		p_portion =  1 - self.probs['p_whole']
		is_whole = bool(self.rng.choice(a = [True, False], p = [self.probs['p_whole'], p_portion]))
		
		# get if integration should be rearranged
		p_not_rearrange = 1 - self.probs['p_rearrange']
		is_rearrange = bool(self.rng.choice(a = [True, False], p = [self.probs['p_rearrange'], p_not_rearrange]))
			
		
		# get if integration should contain deletion
		p_not_delete = 1 - self.probs['p_delete']
		is_delete = bool(self.rng.choice(a = [True, False], p = [self.probs['p_delete'], p_not_delete]))
		
		# get number of fragments - ignored if both isDelete and isRearrange are both False
		
		# must have at least two pieces for a rearrangment, or three for a deletion
		min_split = 0
		if is_rearrange is True:
			min_split = 2
		if is_delete is True:
			min_split = 3
		
		n_fragments = self.poisson_with_minimum_and_maximum(self.probs['lambda_split'], min = min_split)
		if n_fragments is None:
			return
		assert n_fragments >= 0 
		
		# if we're doing a deletion, get the number of fragments to delete
		if is_delete is True:
			n_delete = int(self.rng.choice(range(0, n_fragments - 1)))
		else:
			n_delete = 0
			
		# if we're doing a rearrangement, get the number of swaps to make
		if is_rearrange is True:
			n_swaps = self.poisson_with_minimum_and_maximum(self.probs['lambda_split'], min = 1)
		else:
			n_swaps = 0
				
		# get type of left junction
		p_clean = 1 - self.probs['p_overlap'] - self.probs['p_gap']
		prob_juncs = [self.probs['p_overlap'], self.probs['p_gap'], p_clean]
		junc_types = self.rng.choice(a = ['overlap', 'gap', 'clean'], size = 2, p = prob_juncs)
				
		# get number of bases in left overlap or gap
		n_left_junc = self.poisson_with_minimum_and_maximum(self.probs['lambda_junction'], min = 0)
		n_right_junc = self.poisson_with_minimum_and_maximum(self.probs['lambda_junction'], min = 0)
		if n_left_junc is None or n_right_junc is None:
			return
		
		print(f"trying integration {len(self)}")
		
		# make an integration
		integration = Integration(self.host, 
								  self.virus,
								  rng = self.rng,
								  int_id = len(self),
								  is_whole = is_whole,
								  is_rearrange = is_rearrange,
								  n_swaps = n_swaps,
								  is_delete = is_delete,
								  n_delete = n_delete,
								  n_fragments = n_fragments,
								  junc_types = list(junc_types),
								  n_juncs = (n_left_junc, n_right_junc)
								  )
		
		# append to self if nothing went wrong with this integration
		if integration.chunk.pieces is not None:
			if integration.pos is not None:
				assert integration.id not in [item.id for item in self]
				self.append(integration)
				return True
				
		return False	
		
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
						assert integration.pos is not None
						# write host bases before this integration
						handle.write(str(self.host[chr].seq[position:integration.pos]))
						# write left gap if relevant
						if integration.junc_types[0] == 'gap':
							handle.write(integration.junc_bases[0])
						# write integrated viral bases
						handle.write(integration.chunk.bases)
						# write right gap if relevant
						if integration.junc_types[1] == 'gap':
							handle.write(integration.junc_bases[1])
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
		 - hPos: 0-based position in the original host fasta at which integration occurs
		 	(relative to input fasta)
		 - intStart: position of the ambiguous bases (gap, overlap or clean junction) on left side,
		 	in final output fasta (accounting for any previous integrations)
		 - intStop: position of the ambiguous bases (gap, overlap or clean junction) on left side
		 	in final output fasta (accounting for any previous integrations)
		 - hDeleted: number of bases deleted from host chromosome
		 - virus: name of integrated virus
		 - viral breakpoints: a comma separated list of viral breakpoints which together indicate the integrated bases
			adjacent pairs of breakpoints indicate portions of the virus that have been integrated
		 - vOris: orientation (+ or -) of each portion of the virus that was integrated
		 - juncTypes: types (gap, overlap, clean) of left and right junctions, respectively
		 - juncBases: bases at left and right junctions, respectively
		"""
		# dictionary will keep track of the number of bases previously integrated on each chromosome
		previous_ints = {key:0 for key in self.host.keys()}
		deleted_bases = {key:0 for key in self.host.keys()}
		
		self.sort()
		with open(filename, 'w', newline='') as csvfile:
			intwriter = csv.writer(csvfile, delimiter = '\t')
			intwriter.writerow(['id', 'chr', 'hPos', 'leftStart', 'leftStop',  'rightStart', 'rightStop', 'hDeleted', 'virus', 'vBreakpoints', 'vOris', 'juncTypes', 'juncBases', 'juncLengths', 'whole', 'rearrangement', 'deletion'])
			for i, integration in enumerate(self):
				assert integration.pos is not None
				
				# calculate start and stop position for this integration				
				left_start = integration.pos + previous_ints[integration.chr]
				left_stop = left_start + integration.n_juncs[0]
				
				right_start = left_start + len(integration.chunk.bases)
				right_stop = right_start + integration.n_juncs[1]
				
				# update previous_ints - total integrated bases
				# are the integrated viral bases, and the bases in the gaps/overlaps
				previous_ints[integration.chr] += len(integration.chunk.bases)
				previous_ints[integration.chr] += integration.n_juncs[0]
				previous_ints[integration.chr] += integration.n_juncs[1]

				# format lists into comma-separated strings
				breakpoints = ",".join([str(i) for i in integration.chunk.pieces])
				oris = ','.join(integration.chunk.oris)
				junc_types = ",".join(integration.junc_types)
				junc_bases = ",".join(integration.junc_bases)
				junc_lengths = ",".join([str(i) for i in integration.n_juncs])
				
				
				intwriter.writerow([integration.id,
									integration.chr, 
									integration.pos, 
									left_start, 
									left_stop,
									right_start,
									right_stop, 
									integration.host_deleted,
									integration.chunk.virus,
									breakpoints,
									oris,
									junc_types,
									junc_bases,
									junc_lengths,
									integration.is_whole,
									integration.is_rearrange,
									integration.is_delete])
			
	def __str__(self):
		return f"Viral integrations object with {len(self)} integrations of viral sequences {list(self.virus.keys())} into host chromosomes {list(self.host.keys())}"

	def __repr__(self):
		return f"Object of type Integrations with properties {self}"
		
	
class Integration(dict):
	"""
	Class to store the properties of an individual integration
	"""
	def __init__(self, host, virus, rng, int_id, is_whole = True, is_rearrange = False, n_swaps = 1, is_delete = False, n_delete = 1, n_fragments = 3, junc_types = ('clean', 'clean'), n_juncs = (0, 0), search_length_overlap = 10000):
		"""
		Function to initialise Integration object
		portionType is 'whole' or 'portion' - the part of the virus that has been inserted
		overlapType is two-member tuple of 'clean', 'gap' or 'overlap' - defines the junction at each end of the integration
		
		objects of this class have the following properties:
		self.chr - string: host chromosome on which integration should be located
		self.chunk - ViralChunk object which is integrated
		self.is_whole - boolean indicating if ViralChunk is the whole virus or just a portion
		self.is_rearrange - boolean indicating if ViralChunk is rearranged
		self.is_delete - boolean indicating if ViralChunk is deleted
		self.junc_types - iterable of length 2 indicating the type of overlap at the left and right junctions, respecivley
		junc_types may be 'clean', 'gap' or 'overlap'
		self.n_juncs - iterable of length 2 indicating the number of bases involved in the left and right junctions, respectivley
		
		
		if anything went wrong with initialisation, self.pos is set to None - check this to make sure integration
		is valid
		"""

		
		# check inputs	
		assert isinstance(virus, dict)
		assert is_whole is True or is_whole is False 
		assert is_rearrange is True or is_rearrange is False
		assert is_delete is True or is_delete is False
		assert isinstance(n_swaps, int) and n_swaps >= 0
		assert isinstance(n_delete, int) and n_delete >= 0
		assert is_delete is False or (n_delete < n_fragments - 1)
		assert len(junc_types) == 2
		assert all([i in ['overlap', 'gap', 'clean'] for i in junc_types])
		assert len(n_juncs) == 2
		assert all([i >= 0 for i in n_juncs])
		assert all([isinstance(i, int) for i in n_juncs])
		assert search_length_overlap > 0 and isinstance(search_length_overlap, int)
		
				
		# set parameters that won't be changed
		self.search_length_overlap = search_length_overlap
		self.id = int_id
	
		# get random chromosome on which to do insertion
		self.chr = str(rng.choice(list(host.keys())))
		
		# get viral chunk
		self.chunk = ViralChunk(virus, 
								rng, 
								is_whole = is_whole, 
								deletion = is_delete, 
								n_delete = n_delete, 
								rearrange = is_rearrange, 
								n_swaps = n_swaps, 
								n_fragments = n_fragments)
		assert self.chunk.pieces is not None
		
		# assign properties to this integrations 
		self.is_whole = is_whole
		self.is_rearrange = is_rearrange
		self.is_delete = is_delete
		
		
		# set overlap types
		assert len(junc_types) == 2
		assert all([True if (i in ['clean', 'gap', 'overlap']) else False for i in junc_types])
		self.junc_types = junc_types
		
		# set the number of bases
		self.n_juncs = [n_juncs[0], n_juncs[1]]
		# but overwrite in the case of a clean junction
		if self.junc_types[0] == 'clean':
			self.n_juncs[0] = 0
		if self.junc_types[1] == 'clean':
			self.n_juncs[1] = 0

		assert len(self.n_juncs) == 2
		assert len(self.junc_types) == 2
		
		# number of bases in overlaps must be less than the length of the integrated chunk
		if self.n_overlap_bases() >= len(self.chunk.bases):
			self.pos = None
			return
		
		# set bases belonging to junction
		self.junc_bases = (self.get_junc_bases(rng, 'left'), self.get_junc_bases(rng, 'right'))
				
		# get a position at which to integrate
		if self.get_int_position(host[self.chr].seq, rng) is False:
			self.pos = None
			return
		
		# number of bases deleted from host chromosome
		self.host_deleted = 0 # TODO
		
		# double check for valid chunk

		assert self.chunk.bases == self.chunk.get_bases(virus)
		assert len(self.junc_bases[0]) == self.n_juncs[0]
		assert len(self.junc_bases[1]) == self.n_juncs[1]
		assert all([len(i) == 2 for i in self.chunk.pieces])
		
	def n_overlap_bases(self):
		"""
		Get the total number of bases in overlaps
		"""
		assert len(self.n_juncs) == 2
		assert len(self.junc_types) == 2	
		n = 0
		if self.junc_types[0] == 'overlap':
			n += self.n_juncs[0]
		if self.junc_types[1] == 'overlap':
			n += self.n_juncs[1]
		return n
		
	def get_int_position(self, chr, rng):
		"""
		get a position at which to perform the integration
		the position depends on the overlap type at each side of the overlap
		if both sides are 'clean' or 'gap', position can be random
		'overlap' junctions place constraints on the integration location because the need
		to be placed where there are overlaps
		"""

		
		# if at both ends the junction is either 'clean' or 'gap', just get a random positon
		if all([True if (i in ['clean', 'gap']) else False for i in self.junc_types]):
			self.pos = int(rng.integers(low = 0, high = len(chr)))
			return True

		# if overlap at both ends, look for both overlaps in host chromosome next to each other
		elif self.junc_types[0] == 'overlap' and self.junc_types[1] == 'overlap':
		
			# make string with both overlap bases 
			self.pos = self.find_overlap_bases(self.junc_bases[0] + self.junc_bases[1], chr, rng)
			
			# check for unsuccessful find
			if self.pos == -1:
				self.pos = None
				return False
			
			## need to remove overlapped bases from viral chunk, and adjust chunk start, breakpoints and oris
			self.delete_left_bases(self.n_juncs[0])
			self.delete_right_bases(self.n_juncs[1])
			
			return True
		
		# if one end is an overlap, find those bases in the host chromosome
		# left overlap
		elif self.junc_types[0] == 'overlap':
			
			# find position with overlap at which to do overlap
			self.pos = self.find_overlap_bases(self.junc_bases[0], chr, rng)
			# check for unsuccessful find
			if self.pos == -1:
				self.pos = None
				return False
			
			## need to remove overlapped bases from viral chunk, and adjust chunk start, breakpoints and oris
			self.delete_left_bases(self.n_juncs[0])
			
			# add checks to make sure that chunk is still valid 
			# TODO
			
			return True

		# right overlap
		elif self.junc_types[1] == 'overlap':

			# find position with overlap at which to do overlap
			self.pos = self.find_overlap_bases(self.junc_bases[1], chr, rng)
			# check for unsuccessful find
			if self.pos == -1:
				self.pos = None
				return False
			
			## need to remove overlapped bases from viral chunk, and adjust chunk start, breakpoints and oris
			self.delete_right_bases(self.n_juncs[1])
			
			# add checks to make sure that chunk is still valid 
			# TODO
			
			return True
			
		else:
			raise ValueError(f"junction types {self.junc_types} are not implemented yet")	
			
		return False		
	
	def find_overlap_bases(self, bases, chr, rng):
		"""
		find bases from an overlap in the host chromosome
		"""
		
		# get position around which to search
		start = int(rng.integers(low = 0, high = len(chr)))
		
		# get start and stop positions for bases in host chromosome to search for overlaps
		stop = start + self.search_length_overlap
		
		# make sure that we aren't searching after the end of the chromosome
		if stop > len(chr):
			stop = len(chr)
			
		found = re.search(bases, str(chr[start:stop]), re.IGNORECASE)
			
		if not found:
		# check for unsuccessful find
			return -1
			
		return found.span()[0] + start
			
	def delete_left_bases(self, n):
		"""
		delete bases on the left after adding an overlap - need to adjust chunk bases, oris and breakpoints and start
		"""
		
		# check we're not trying to delete more bases than there are in the chunk
		assert n < len(self.chunk.bases)
		
		# adjust bases
		self.chunk.bases = self.chunk.bases[n:]
		
		# adjust breakpoints - delete bases one by one from breakpoints
		# until we've deleted enough baes
		deleted_bases = 0
		to_delete = []
		i = 0
		while deleted_bases < n:
			# if we're left with a piece of length 0, flag this piece for deletion
			if self.chunk.pieces[i][0] == self.chunk.pieces[i][1]:
				to_delete.append(i)
				i += 1
			# if this piece is a forward piece
			if self.chunk.oris[i] == '+':
				# detele one base
				self.chunk.pieces[i][0] += 1
				deleted_bases += 1
			# if this piece is a reverse piece we need to remove from the end
			# because self.chunk.bases has already taken orientations into account
			elif self.chunk.oris[i] == "-":
				self.chunk.pieces[i][1] -= 1
				deleted_bases += 1
			else:
				print(f"unrecgonised orientation {self.chunk.oris[i]} in chunk {vars(self.chunk)}")
				self.pos = None
				return
				
		# remove chunks that we want to dele
		self.chunk.pieces = [self.chunk.pieces[i] for i in range(len(self.chunk.pieces)) if (i not in to_delete)]
		self.chunk.oris = [self.chunk.oris[i] for i in range(len(self.chunk.oris)) if (i not in to_delete)]
		
		#adjust start and stop
		breakpoints = [piece[0] for piece in self.chunk.pieces]
		breakpoints += [piece[1] for piece in self.chunk.pieces]
		breakpoints.sort()
		self.chunk.start = breakpoints[0]
		self.chunk.stop = breakpoints[-1]
		
	def delete_right_bases(self, n):
		"""
		delete bases on the left after adding an overlap - need to adjust chunk bases, oris and breakpoints and stop
		"""
		
		# check we're not trying to delete more bases than there are in the chunk
		assert n < len(self.chunk.bases)
		
		# adjust stop
		self.chunk.stop -= n
		
		# adjust bases
		self.chunk.bases = self.chunk.bases[:-n]
		
		# adjust breakpoints 		
		deleted_bases = 0
		to_delete = []
		i = 0
		while deleted_bases < n:
			# if we're left with a piece of length 0, flag this piece for deletion
			if self.chunk.pieces[i][0] == self.chunk.pieces[i][1]:
				to_delete.append(i)
				i += 1
			# if this piece is a forward piece
			if self.chunk.oris[i] == "+":
				# delete one base
				self.chunk.pieces[i][1] -= 1
				deleted_bases += 1
			elif self.chunk.oris[i] == "-":
				# delete one base
				self.chunk.pieces[i][0] += 1
				deleted_bases += 1
			else:
				print(f"unrecgonised orientation {self.chunk.oris[i]} in chunk {vars(self.chunk)} ")
				self.pos = None
				return
				 
		# remove chunks that we want to delete
		self.chunk.pieces = [self.chunk.pieces[i] for i in range(len(self.chunk.pieces)) if (i not in to_delete)]
		self.chunk.oris = [self.chunk.oris[i]  for i in range(len(self.chunk.oris)) if (i not in to_delete)]
	
		#adjust start and stop
		breakpoints = [piece[0] for piece in self.chunk.pieces]
		breakpoints += [piece[1] for piece in self.chunk.pieces]
		breakpoints.sort()
		self.chunk.start = breakpoints[0]
		self.chunk.stop = breakpoints[-1]
			
	def get_junc_bases(self, rng, side):
		"""
		Get the bases at the left or right junction, depending on whether the junction is
		a gap, overlap, or clean junction
		"""
		assert side in ['left', 'right']
		assert len(self.junc_types) == 2
		assert len(self.n_juncs) == 2
		assert self.junc_types[0] in ['gap', 'overlap', 'clean']
		assert self.junc_types[1] in ['gap', 'overlap', 'clean']

		
		if side == 'left':
			n_bases = self.n_juncs[0]
			# no bases in a clean junction
			if self.junc_types[0] == 'clean':
				return ""
			# random bases in a gap
			elif self.junc_types[0] == 'gap':
				return self.get_n_random_bases(rng, n_bases)
			# first n bases of viral chunk in an overlap
			elif self.junc_types[0] == 'overlap':
				return self.chunk.bases[:n_bases]
			else:
				raise ValueError(f"unrecgonised type: {self.junc_types[0]}")
		elif side == 'right':
			n_bases = self.n_juncs[1]
			if self.junc_types[1] == "clean":
				return ""
			elif self.junc_types[1] == 'gap':
				return self.get_n_random_bases(rng, n_bases)
			elif self.junc_types[1] == 'overlap':
				return self.chunk.bases[-n_bases:]
			else:
				raise ValueError(f"unrecgonised type: {self.junc_types[1]}")
		else:
			raise ValueError(f"unrecgonised side: {side}")	
	
	def get_n_random_bases(self, rng, n_bases):
		"""
		get a string composed of n random nucleotides
		"""	
		return "".join(rng.choice(['A', 'T', 'G', 'C'], n_bases))
		
	def __str__(self):
		return f"Viral integration into host chromosome {self.chr}'"
		
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

	
class ViralChunk(dict):
	"""
	Class to store properties of an integrated chunk of virus
	"""
	
	def __init__(self, virus,  rng, is_whole = False, deletion=True, n_delete = 1, rearrange=True, n_swaps = 1, n_fragments = 3):
		"""
		function to get a chunk of a virus
		virus is the dictionary of seqrecords imported using biopython
		is_whole specifies if we want the whole virus (if True), or just a portion (if False)
		
		a viral chunk can be subdivided into smaller pieces if we want
		a deletion or a rearrangement.  if deletion is true, the chunk
		will be subdivided into at least three pieces and one of the middle ones
		will be deleted
		if rearrange is True, the chunk will be subdivided into at least two pieces
		and they will be randomly shuffled

		the bases attribute of a ViralChunk consist of only the bases that are unique to the virus. 
		So in the case of an Integration of a ViralChunk with a 'overlap' type junction,
		the bases, breakpoints and oris attributes are re-assigned to remove the overlapped bases
		
		"""
		# check inputs
		assert isinstance(virus, dict)
		
		# get a random virus to integrate
		self.virus = str(rng.choice(list(virus.keys())))
		
		if is_whole is True:
			self.start = 0
			self.stop = len(virus[self.virus])
		elif is_whole is False:
			self.start = int(rng.integers(low = 0, high = len(virus[self.virus].seq) - 1))
			self.stop = int(rng.integers(low = self.start + 1, high = len(virus[self.virus].seq)))
		else:
			raise ValueError('is_whole must be either True or False')
			
		# breakpoints are the start and stop coordinates of pieces of the virus that have been 
		# integrated
		
		# set breakpoints
		self.pieces = [[self.start, self.stop]]
		
		self.oris = str(rng.choice(('+', '-')))
		
		# do a deletion if applicable
		if deletion is True:
			self.__delete(rng, n_fragments, n_delete)
			# if something went wrong, breakpoints will be None
			if self.pieces is None:
				return
		elif deletion is False:
			pass
		else:
			raise valueError("deletion must be either True or False")
				
		if rearrange is True:
			self.__rearrange(rng, n_fragments, n_swaps)
			# if something went wrong, breakpoints will be None
			if self.pieces is None:
				return
		elif rearrange is False:
			pass
		else:
			raise valueError("rearrange must be either True or False")
		
		# get bases
		self.bases = self.get_bases(virus)
				
	def get_bases(self, virus):
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
		
	def __split_into_pieces(self, rng, n_fragments):
		"""
		get random, unique breakpoints to divide a viral chunk into pieces
		there must be at least min_breakpoints, which results in min_breakpoints + 1 pieces 
		this is a list of coordinates, not tuples (unlike self.pieces)
		"""
		# shouldn't do this if this chunk has already been split
		assert len(self.pieces) == 1
		assert len(self.oris) == 1
		
		# check that we're not trying to divide a chunk into more pieces than there are bases
		if n_fragments >= self.stop - self.start:
			self.pieces = None
			return
		
		# get the number of pieces to divide into
		num_breakpoints = n_fragments - 1
				
		# get random breakpoints from within this chunk
		breakpoints = rng.choice(range(self.start + 1, self.stop - 1), size = num_breakpoints, replace = False)
				
		# set self.pieces
		breakpoints = [self.start] + sorted(breakpoints) + [self.stop]
		self.pieces = [[int(breakpoints[i]), int(breakpoints[i+1])] for i in range(len(breakpoints) - 1)]
		
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
	
	def __delete(self, rng, n_fragments, n_delete):
		"""
		Divide a viral chunk up into multiple pieces
		and remove one of those pieces
		"""
		# deletions are always performed first, so the chunk should not have been split yet
		assert len(self.pieces) == 1
		
		assert n_fragments - n_delete >= 2 # want to have at least two pieces left

		# split chunk into at n_fragments pieces
		self.__split_into_pieces(rng, n_fragments)
		if self.pieces is None:
			return
		assert len(self.pieces) == n_fragments
		
		# decide which portions to delete
		i_delete = rng.choice(range(1, len(self.pieces) - 1), n_delete, replace=False)
		
		# do deletion
		self.pieces = [piece for i, piece in enumerate(self.pieces) if i not in i_delete]
		self.oris = [ori for i, ori in enumerate(self.oris) if i not in i_delete]	
		
		assert len(self.pieces) == n_fragments - n_delete
		
	def __rearrange(self, rng, n_fragments, n_swaps):
		"""
		Divide a viral chunk up into multiple pieces
		and randomise their order and orientiations
		"""
		# split the chunk if it hasn't already been split
		if len(self.pieces) == 1:
			# split chunk into at least three pieces
			self.__split_into_pieces(rng, n_fragments)
			if self.pieces is None:
				return
			assert len(self.pieces) == n_fragments
		else:
			assert len(self.pieces) > 1
		
		# if we only have two pieces, we should only do one swap
		# so that we don't end up back with the same fragment
		# there are other ways to end up with the same fragment after swaps
		# but don't worry about them for now - TODO
		if len(self.pieces) == 2:
			n_swaps = 1
		
		for i in range(n_swaps):
			# pick a point about which to swap
			if 1 == len(self.pieces) - 1:
				i_swap = 1
			else:
				i_swap = rng.choice(range(1, len(self.pieces) - 1))
				
			# swap everything to the left of this position with everything on the right
			self.pieces = self.pieces[i_swap:] + self.pieces[:i_swap]
			
			# 50 % chance of swapping the orientations of all the pieces for each side
			if bool(rng.choice((True, False))) is True:
				self.__swap_orientations(self, i_swap, side = 'left')
			if bool(rng.choice((True, False))) is True:
				self.__swap_orientations(self, i_swap, side = 'right')

	
	def __str__(self):
		return f"Viral chunk of virus {self.virus} ({self.start}, {self.stop}) and orientations {self.oris}"
		
	def __repr__(self):
		return f"Object of type ViralChunk with properties {self}"

if __name__ == "__main__":
	main(argv[1:])
