#### simulates viral integrations ####

#### written by Suzanne Scott (suzanne.scott@csiro.au) and Susie Grigson (susie.grigson@flinders.edu) ####

# types of possible integrations:
## whole viral genome (possibility for reverse orientation)
## portion of viral genome (possibility for reverse orientation)
## n sequential portions of viral genome, rearranged (with possibility that some are reversed)
## n non-sequential (with a gap in between) portions of viral genome, rearranged (with possibility that some are reversed) (ie deletion)
## portion of viral genome divided into n sequential portions of viral genome, rearranged 
## poriton of viral genome divided into n non-sequential portions of viral genome and rearranged 

# at each integration type, overlaps possible are gap, overlap, none
## if gap, insert small number of random bases between host and virus
## if overlap, take small (<=10) number of bases and look for homology, do integration there
## if none, there is a 'clean' junction between virus and host

# reports parameters of integrations:
	#location in host genome (chr, start, stop)
	#part of viral genome inserted (virus, start, stop)
	#integration type (whole, portion, rearrangement, etc)
	#overlaps/gaps at each junction 
# reports all integration sites relative to the host genome, independent of other integrations
# intergrations are stored internally though the Integration class

###import libraries
from Bio import SeqIO
from Bio.Alphabet.IUPAC import unambiguous_dna
from Bio.Seq import Seq 
from Bio.SeqRecord import SeqRecord
from scipy.stats import poisson
import pandas as pd
import argparse
from sys import argv
from os import path
import numpy as np
import scipy
import pdb 
import csv
import re
import pprint
import time

###
max_attempts = 200000 #maximum number of times to try to place an integration site 

default_ints = 5
default_epi = 0

default_p_whole = 0.3
default_p_rearrange = 0.05
default_p_delete = 0.05
default_lambda_split = 1.5
default_p_overlap = 0.3
default_p_gap = 0.3
default_lambda_junction = 3
default_p_host_del = 0.0
default_lambda_host_del = 1000
default_min_sep = 500

search_length_overlap = 10000 # number of bases to search either side of randomly generated position for making overlaps

lambda_junc_percentile = 0.99

fasta_extensions = [".fa", ".fna", ".fasta"]

pp = pprint.PrettyPrinter(indent=4)

def main(argv):
	#get arguments
	parser = argparse.ArgumentParser(description='simulate viral integrations')
	parser.add_argument('--host', help='host fasta file', required = True, type=str)
	parser.add_argument('--simug_snp_indel', help='output file from simuG to map location of snps and indels')
	parser.add_argument('--virus', help = 'virus fasta file', required = True, type=str)
	parser.add_argument('--ints', help = 'output fasta file', type=str, default="integrations.fa")
	parser.add_argument('--int_info', help = 'output tsv with info about viral integrations', type=str, default="integrations.tsv")
	parser.add_argument('--int_num', help = f"number of integrations to be carried out [{default_ints}]", type=int, default=default_ints)
	parser.add_argument('--p_whole', help = f"probability of a virus being integrated completely [{default_p_whole}]", type = float, default=default_p_whole) 
	parser.add_argument('--p_rearrange', help=f"probability of an integrated piece of virus being rearranged [{default_p_rearrange}]", type=float, default=default_p_rearrange) 
	parser.add_argument('--p_delete', help=f"probability of an integrated piece of virus containing a deletion [{default_p_delete}]", type=float, default=default_p_delete) 
	parser.add_argument('--lambda_split', help = f"mean of poisson distriubtion for number of pieces to split an integrated viral chunk into when rearranging or deleting [{default_lambda_split}]", type=float, default=default_lambda_split)
	parser.add_argument('--p_overlap', help=f"probability of a junction to be overlapped (common bases between virus and host) [{default_p_overlap}]", type=float, default=default_p_overlap) 
	parser.add_argument('--p_gap', help=f"probability of a junction to have a gap [{default_p_gap}]", type=float, default=default_p_gap) 
	parser.add_argument('--lambda_junction', help = f"mean of poisson distribution of number of bases in a gap/overlap [{default_lambda_junction}]", type=float, default=default_lambda_junction)
	parser.add_argument('--p_host_deletion', help=f"probability of a host deletion at the integation site [{default_p_host_del}]", type=float, default=default_p_host_del) 	
	parser.add_argument('--lambda_host_deletion', help = f"mean of poisson distribution of number of bases deleted from the host [{default_lambda_host_del}]", type=float, default=default_lambda_host_del)	
	parser.add_argument('--epi_num', help = f"number of episomes to added to output fasta [{default_ints}]", type=int, default=default_ints)
	parser.add_argument('--epi_info', help = 'output tsv with info about episomes', type=str, default="episomes.tsv")	
	parser.add_argument('--seed', help = 'seed for random number generator', default=12345, type=int)
	parser.add_argument('--verbose', help = 'display extra output for debugging', action="store_true")
	parser.add_argument('--model-check', help = 'check integration model every new integration', action="store_true")
	parser.add_argument('--min-sep', help='minimum separation for integrations', type=int, default=default_min_sep)
	parser.add_argument('--min-len', help='minimum length of integrated viral fragments', type=int, default=None)

	args = parser.parse_args()
	
	#### generate integrations ####		
	
	# generate dictionary of probabilities and lambda/means as input for Integrations class
	probs = {'p_whole' 				: args.p_whole,
			'p_rearrange' 			: args.p_rearrange,
			'p_delete'				: args.p_delete,
			'lambda_split'			: args.lambda_split,
			'p_overlap' 			: args.p_overlap,
			'p_gap' 				: args.p_gap,
			'lambda_junction' 		: args.lambda_junction,
			'p_host_del'			: args.p_host_deletion,
			'lambda_host_del'		: args.lambda_host_deletion}
			
	# initialise Events
	if args.verbose is True:
		print("initialising a new Events object")

	seqs = Events(args.host, args.virus, seed=args.seed, verbose = True, 
								min_len = args.min_len, simug_snp_indel = args.simug_snp_indel)
	
	# add integrations
	seqs.add_integrations(probs, args.int_num, max_attempts, 
												model_check = args.model_check, min_sep = args.min_sep)
	seqs.add_episomes(probs, args.epi_num, max_attempts)
	
	
	# save outputs
	seqs.save_fasta(args.ints)
	seqs.save_integrations_info(args.int_info)
	seqs.save_episomes_info(args.epi_info)
	
class Events(dict):
	"""
	base class for this script - stores two kinds of events: Integrations and Episomes
	Integrations is a list-like class of Integration objects, which contain information about pieces of 
	virus that have been integrated and their junctions with the surrounding chromomsomal sequence
	Episomes is a list-like class of ViralChunk objects, which contain information about pieces of virus that
	are episomal (present in sequence data but not integrated into the host chromosomes)
	"""
	
	def __init__(self, host_fasta_path, virus_fasta_path, 
								fasta_extensions = ['.fa', '.fna', '.fasta'] , 
								seed = 12345, verbose = False, min_len = None,
								simug_snp_indel = None, simug_cnv = None, 
								simug_inversion = None, simug_tranlocation = None):
		"""
		initialise events class by importing a host and viral fasta, and setting up the random number generator
		"""
		
		# expectations for inputs
		assert isinstance(fasta_extensions, list)
		assert isinstance(seed, int)
		assert isinstance(verbose, bool)
		assert isinstance(min_len, int) or min_len is None
		if min_len is not None:
			assert min_len > 0
		
		self.verbose = verbose
		self.min_len = min_len
		
		# check and import fasta files
		if self.verbose is True:
			print("importing host fasta", flush=True)
	
		# read host fasta - use index which doesn't load sequences into memory because host is large genome
		if self.checkFastaExists(host_fasta_path, fasta_extensions):
			self.host = SeqIO.to_dict(SeqIO.parse(host_fasta_path, 'fasta', alphabet=unambiguous_dna))
		else:
			raise OSError("Could not open host fasta")
		
		if self.verbose is True:
			print(f"imported host fasta with {len(self.host)} chromosomes:", flush=True)
			host_chrs = {key:len(self.host[key]) for key in self.host.keys()}
			for key, length in host_chrs.items():
				print(f"\thost chromosome '{key}' with length {length}", flush=True)
	
		# read virus fasta -  make dictionary in memory of sequences
		if self.checkFastaExists(virus_fasta_path, fasta_extensions):
			self.virus = SeqIO.to_dict(SeqIO.parse(virus_fasta_path, 'fasta', alphabet=unambiguous_dna))
			# convert all viral sequences to upper case to make it easier to check output
			for this_virus in self.virus.values():
				this_virus.seq = this_virus.seq.upper()
		else:
			raise OSError("Could not open virus fasta")
			
		# check for SNP and indel map file from simuG - need to account for indels in output
		self.indel = [row for row in self.read_simug_file(simug_snp_indel) if row['variant_type'] == 'INDEL']
		
		# other types of structural variants
		self.cnv = self.read_simug_file(simug_cnv)
		self.inversion = self.read_simug_file(simug_inversion)
		self.tranlocation = self.read_simug_file(simug_tranlocation)

		# check that minimum length is greater than the length of all the viral references
		if self.min_len is not None:
			if not all([len(virus) >= self.min_len for virus in self.virus.values()]):
				raise ValueError(f"specified minimum length is more than the length of one of the input viruses")
			if any([len(virus) == self.min_len for virus in self.virus.values()]):
				print(f"warning: minimum length is equal to the length of one or more references - integrations involving these references will all be whole, regardless of p_whole")


		if self.verbose is True:
			print(f"imported virus fasta with {len(self.virus)} sequences:", flush=True)
			virus_seqs = {key:len(self.virus[key]) for key in self.virus.keys()}
			for key, length in virus_seqs.items():
				print(f"\tviral sequence '{key}' with length {length}", flush=True)
			
		# instantiate random number generator
		self.rng = np.random.default_rng(seed)
		
	def read_simug_file(self, filename):
		"""
		open simuG output file and read contents into memory
		"""
		if filename is None:
			return None
			
		assert path.isfile(filename)
		lines = []
		with open(filename, 'r') as infile:
			reader = csv.DictReader(infile, delimiter =  '\t')
			for row in reader:
				lines.append(row)
				
		return lines
		 
	
	def add_integrations(self, probs, int_num, max_attempts = 50, model_check = False, min_sep = 1):
		"""
		Add an Integrations object with int_num integrations, with types specified by probs,
		to self
		"""
		
		assert isinstance(max_attempts, int)
		assert isinstance(int_num, int)
		assert isinstance(min_sep, int)
		assert max_attempts > 0
		assert int_num >= 0
		assert min_sep > 0
		
		self.min_sep = min_sep
		self.max_attempts = max_attempts
		self.model_check = model_check
		
		# can only add integrations once
		if 'ints' in self:
			raise ValueError("integrations have already been added to this Events object") # is there a better error type for this?
			
		# check that the number of requested integrations will fit in the host, given the requested minimum separation
		# rule of thumb is that we allow 4*min_sep for each integration
		total_host_length = sum([len(seq) for seq in self.host.values()])
		if self.min_sep * 4 * int_num > total_host_length:
			raise ValueError("The requested number of integrations, with the specified minimum separation, are not likely to fit into the specified host.  Either decrease the number of integrations or the minimum separation.")
			
		# we require that the minimum length of integrations is longer than the integrated virus
		# so check the value of lambda_junction relative to min_len and the length of the shortest viral reference
		# require that both are greater than the 99th percentile of the poisson distribution defined by lambda_junction
		self.check_junction_length(probs)
			
		# instantiate Integrations object
		self.ints = Integrations(self.host, self.virus, probs, self.rng, self.max_attempts, self.model_check, self.min_sep, self.min_len, self.indel)
		
		# add int_num integrations
		if self.verbose is True:
			print(f"performing {int_num} integrations", flush=True)
		counter = 0
		while len(self.ints) < int_num:
			t0 = time.time()
			if self.ints.add_integration() is False:
				counter += 1
			# check for too many attempts
			if counter > max_attempts:
				raise ValueError('too many failed attempts to add integrations')
			t1 = time.time()
			print(f"added integration {len(self.ints)} in {(t1-t0):.2f}s", flush=True)
			print()
		# if we had fewer than 50% of our attempts left
		if (counter / max_attempts) > 0.5:
			print(f"warning: there were {counter} failed integrations", flush=True)
				
	def add_episomes(self, probs, epi_num, max_attepmts = 50):
		"""
		Add an Integrations object with int_num integrations, with types specified by probs,
		to self
		"""
		
		assert isinstance(max_attempts, int)
		assert isinstance(epi_num, int)
		assert max_attempts > 0
		assert epi_num >= 0
		
		# can only add episomes once
		if 'epis' in self:
			raise ValueError("episomes have already been added to this Events object")
			
		# we require that the minimum length of integrations is longer than the integrated virus
		# so check the value of lambda_junction relative to min_len and the length of the shortest viral reference
		# require that both are greater than the 99th percentile of the poisson distribution defined by lambda_junction
		self.check_junction_length(probs)
			
		# instantiate Episomes object
		self.epis = Episomes(self.virus, self.rng, probs, max_attempts, self.min_len)
		
		# add epi_num episomes
		if self.verbose is True:
			print(f"adding {epi_num} episomes", flush=True)
		counter = 0
		while len(self.epis) < epi_num:
			if self.epis.add_episome() is False:
				counter += 1
			if counter > max_attempts:
				raise ValueError('too many failed attempts to add episomes')

			
	def checkFastaExists(self, file, fasta_extensions):
		#check file exists
		exists = path.isfile(file)
		if not(exists):
			return False
		#check extension
		prefix = path.splitext(file)[-1]
		if prefix:
			if not prefix in fasta_extensions:
				return False
		return True		
		
	def check_junction_length(self, probs):
		"""
		we require that the minimum length of integrations is longer than the integrated virus
		so check the value of lambda_junction relative to min_len and the length of the shortest viral reference
		require that both are greater than the 99th percentile of the poisson distribution defined by lambda_junction
		"""
		if probs['p_gap'] + probs['p_overlap'] == 0:
			return
		thresh = poisson.ppf(lambda_junc_percentile, probs['lambda_junction'])
		if self.min_len is not None:
			if thresh * 2  > self.min_len:
				raise ValueError(
					"There is likely to be a lot of clashes between the length of the left and right junctions, and the \
					length of the integrations.  Set a shorter lambda_jucntion or a longer minimum length"
					)
		if thresh * 2  > min([len(virus) for virus in  self.virus.values()]):
			raise ValueError(
				"There is likely to be a lot of clashes in which the length of the left and right junctions, and the \
				length of the integrations.  Set a shorter lambda_jucntion or use longer viral references."
				)
		
	def save_fasta(self, file):
		"""
		save output sequences to file
		"""

		if 'ints' in vars(self):
			assert len(self.ints) >= 0
			self.ints._Integrations__save_fasta(file, append = False)
		
		if 'epis' in vars(self):
			assert len(self.epis) >= 0
			self.epis._Episomes__save_fasta(file, append = True)
			
		if ('ints' not in vars(self)) and ('epis' not in vars(self)):
			print("warning: no integrations or episomes have been added")
			
		if self.verbose is True:
			print(f"saved fasta with integrations and episomes to {file}")
			
	def save_integrations_info(self, file):
		"""
		save info about integrations to file
		"""
		assert len(self.ints) >= 0
		
		self.ints._Integrations__save_info(file)
		
		if self.verbose is True:
			print(f"saved information about integrations to {file}")
		
	def save_episomes_info(self, file):
		"""
		save info about episomes to file
		"""
		assert 'epis' in vars(self)
		assert len(self.epis) >= 0
		
		self.epis._Episomes__save_info(file)
		
		if self.verbose is True:
			print(f"saved information about episomes to {file}", flush=True)
		
class Integrations(list):
	"""
	Class to store all integrations for a given host and virus
	This class stores the sequences for the host and virus (dictionaries of biopython seqRecord objects)
	And probabilities of each type of viral integration - whole/portion, rearrangement/deletion, gap/overlap, 
	means of poisson distributions, host deletions and their mean lengths, etc
	
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
	
	p_host_del - probability that there will be a deletion from the host at integration site
	lambda_host_del - mean of poisson distribution of number of bases deleted from host genome at integration site

	"""
	def __init__(self, host, virus, probs, rng, max_attempts=50, model_check=False, min_sep=1, min_len=None, indel=None):
		"""
		Function 
		to initialise Integrations object
		"""
		
		# checks for inputs
		assert isinstance(virus, dict)
		assert isinstance(probs, dict)
		assert isinstance(max_attempts, int)
		assert max_attempts > 0
		assert isinstance(model_check, bool)
		assert isinstance(min_sep, int)
		assert min_sep > 0
		assert isinstance(min_len, int) or min_len is None
		if min_len is not None:
			assert min_len > 0
		assert isinstance(indel, list) or indel is None
		
		# assign properties common to all integrations
		self.host = host
		self.virus = virus
		self.probs = probs
		self.rng = rng
		self.max_attempts = max_attempts
		self.model_check = model_check
		self.min_sep = min_sep
		self.min_len = min_len
		self.indel = indel
		
		# default values for probs
		self.set_probs_default('p_whole', default_p_whole)
		self.set_probs_default('p_rearrange', default_p_rearrange)
		self.set_probs_default('p_delete', default_p_delete)
		self.set_probs_default('lambda_split', default_lambda_split)
		self.set_probs_default('p_overlap', default_p_overlap)
		self.set_probs_default('p_gap', default_p_gap)
		self.set_probs_default('lambda_junction', default_lambda_junction)	
		self.set_probs_default('p_host_del', default_p_host_del)
		self.set_probs_default('lambda_host_del', default_lambda_host_del)	
		
		# checks on assigned variables
		assert 'p_whole' in probs
		assert 'p_rearrange' in probs
		assert 'p_delete' in probs	
		assert 'lambda_split' in probs	
		assert 'p_overlap' in probs
		assert 'p_gap' in probs
		assert 'lambda_junction' in probs
		assert 'p_host_del' in probs
		assert 'lambda_host_del' in probs
					
		# check that probabilities are between 0 and 1
		self.check_prob(self.probs['p_whole'], 'p_whole')
		#self.check_prob(self.probs['p_rearrange'] + self.probs['p_delete'], 'the sum of p_rearrange and p_delete')	
		self.check_prob(self.probs['p_rearrange'], 'p_rearrange')	
		self.check_prob(self.probs['p_delete'], 'p_delete')	
		self.check_prob(self.probs['p_overlap'] + self.probs['p_gap'], 'the sum of p_overlap and p_gap')
		
		self.check_prob(self.probs['p_host_del'], 'p_host_deletion')
			
		# check that lambda values are positive floats	
		self.check_float_or_int_positive(self.probs['lambda_split'], 'lambda_split')
		self.check_float_or_int_positive(self.probs['lambda_junction'], 'lambda_junction')
		self.check_float_or_int_positive(self.probs['lambda_host_del'], 'lambda_host_deletion')
		
		# initialize model of integrations
		# each chromosome is an entry in the model dict
		# each chromosome is composed of a list of sequences, each composed of a dictionary
		# each sequence is a dictionary, specifying the 'origin' of the bases - 'host', 'virus', 'ambig'
		# as well as the start and stop and orientation of that sequence (in a list)
		# this is to be updated every time we do an integration

		self.model = {chr : [{'origin': 'host', 'coords':(0, len(seq)), "ori" : "+", 'seq_name': chr}] for chr, seq in host.items()}
			

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
				print(f"warning: too many attempts to get value with minimum {min} and maximum {max} from poisson distribution with mean {lambda_poisson}", flush=True)
				
				return None
					
	def add_integration(self):
		"""
		Add an integration by appending an Integration object to self.  
		"""
			
		# call functions that randomly set properties of this integrations
		counter = 0
		while True:
			chunk_props = self.set_chunk_properties()
			if len(chunk_props) != 0:
				break
			counter += 1
			if counter > self.max_attempts:
				raise ValueError("too many attempts to set chunk properties")
		assert len(chunk_props) == 5
			
		counter = 0
		while True:
			junc_props = self.set_junc_properties()
			if len(junc_props) != 0:
				break
			counter += 1
			if counter > self.max_attempts:
				raise ValueError("too many attempts to set junction properties")
		
		assert len(junc_props) == 4

		# make an integration
		integration = Integration(self.host, 
								  self.virus,
								  model = self.model,
								  rng = self.rng,
								  int_id = len(self),
								  chunk_props = chunk_props,
								  junc_props = junc_props,
								  min_sep = self.min_sep
								  )
		
		# append to self if nothing went wrong with this integration
		if integration.chunk.pieces is not None:
			if integration.pos is not None:
				assert integration.id not in [item.id for item in self]
				
				self.append(integration)
				
				self.__update_model(integration)
				
				return True
				
		return False
		
	def set_junc_properties(self):
		"""
		randomly set the properties of the junctions (self.junc_props) for this Integration
		dict with the following properties:
			- junc_types = iterable of length 2 specifying type of left and right junctions (one of gap, overlap, clean)
			- n_junc = iterable of length 2 specifying length of left and right junctions
			- host_del = boolean specifying if there should be a deletion from the host at the integration site
			- host_del_len = integer specifying the number of bases to be deleted from the host at the integration site
		"""
		
		junc_props = {}
		
		# get type of left junction
		p_clean = 1 - self.probs['p_overlap'] - self.probs['p_gap']
		prob_juncs = [self.probs['p_overlap'], self.probs['p_gap'], p_clean]
		junc_props['junc_types'] = self.rng.choice(a = ['overlap', 'gap', 'clean'], size = 2, p = prob_juncs)
				
		# get number of bases in left junction
		if junc_props['junc_types'][0] == 'clean':
			n_left_junc = 0
		elif junc_props['junc_types'][0] in ['gap', 'overlap']:
			n_left_junc = self.poisson_with_minimum_and_maximum(self.probs['lambda_junction'], min = 1)
			if n_left_junc is None:
				return {}
		else:
			return {}
			
		# get number of bases in right junction
		if junc_props['junc_types'][1] == 'clean':
			n_right_junc = 0
		elif junc_props['junc_types'][1] in ['gap', 'overlap']:
			n_right_junc = self.poisson_with_minimum_and_maximum(self.probs['lambda_junction'], min = 1)
			if n_right_junc is None:
				return {}
		else:
			return {}
		
		junc_props['n_junc'] = (n_left_junc, n_right_junc)
		
		# check that if we have a clean junction, it's length is 0 
		assert not(junc_props['junc_types'][0] == 'clean') or junc_props['n_junc'][0] == 0
		assert not(junc_props['junc_types'][1] == 'clean') or junc_props['n_junc'][1] == 0 
		
		# check that if we don't have a clean junction, it's length is greater than zero
		assert not(junc_props['junc_types'][0] != 'clean') or junc_props['n_junc'][0] > 0
		assert not(junc_props['junc_types'][1] != 'clean') or junc_props['n_junc'][1] > 0 
		
		# decide if this integration should have a deletion from the host
		host_deletion = self.rng.choice(a = [True, False], p = (self.probs['p_host_del'], 1 - self.probs['p_host_del']))
		junc_props['host_del'] = bool(host_deletion)
		
		# if we're doing a host deletion, get number of bases to be deleted
		if junc_props['host_del'] is True:
			junc_props['host_del_len'] = self.poisson_with_minimum_and_maximum(self.probs['lambda_host_del'], min = 1)
			if junc_props['host_del_len'] is None:
				return {}
		elif junc_props['host_del'] is False:
			junc_props['host_del_len'] = 0
		else:
			return {}
			
		
		# check that 
		return junc_props
		
	def set_chunk_properties(self):
		"""
		randomly set the properties of the viral chunk for this Integration
		returns dict with the following properties:
			- is_whole: boolean specifying if the ViralChunk is whole (if false, chunk is just a portion)
			- n_fragments: number of fragments into which ViralChunk should be split
			- n_delete: number of fragments to delete from ViralChunk (should always leave at least two fragments after deletion)
			- n_swaps: number of swaps to make when rearranging ViralChunk
			- min_len: the minimum length of the viral chunk - optional (if not present, for integration will be set to
							the number of overlap bases + 1, for episome will be set to 1
		"""
		
		chunk_props = {}
		
		# get if integration should be whole or portion
		p_portion =  1 - self.probs['p_whole']
		chunk_props['is_whole'] = bool(self.rng.choice(a = [True, False], p = [self.probs['p_whole'], p_portion]))
		
		# get if integration should be rearranged
		p_not_rearrange = 1 - self.probs['p_rearrange']
		is_rearrange = bool(self.rng.choice(a = [True, False], p = [self.probs['p_rearrange'], p_not_rearrange]))
			
		
		# get if integration should contain deletion
		p_not_delete = 1 - self.probs['p_delete']
		is_delete = bool(self.rng.choice(a = [True, False], p = [self.probs['p_delete'], p_not_delete]))
		
		# get number of fragments - ignored if both isDelete and isRearrange are both False
		# must have at least two pieces for a rearrangment, or three for a deletion
		min_split = 1
		if is_rearrange is True:
			min_split = 2
		if is_delete is True:
			min_split = 3
		
		# set the number of fragments for the chunk
		if is_delete is False and is_rearrange is False:
			chunk_props['n_fragments'] = 1
		else:
			chunk_props['n_fragments'] = self.poisson_with_minimum_and_maximum(self.probs['lambda_split'], min = min_split)
		if chunk_props['n_fragments'] is None:
			return {}
		assert chunk_props['n_fragments'] > 0 
		
		# if we're doing a deletion, get the number of fragments to delete
		if is_delete is True:
			chunk_props['n_delete'] = int(self.rng.choice(range(0, chunk_props['n_fragments'] - 1)))
		else:
			chunk_props['n_delete'] = 0
			
		# if we're doing a rearrangement, get the number of swaps to make
		if is_rearrange is True:
			chunk_props['n_swaps'] = self.poisson_with_minimum_and_maximum(self.probs['lambda_split'], min = 1)
			if chunk_props['n_swaps'] is None:
				return {}
		else:
			chunk_props['n_swaps'] = 0
			
		# set minimum length of chunk
		chunk_props['min_len'] = self.min_len
			
		return chunk_props
		
	def __update_model(self, integration):
		"""
		update self.model for a new integration 
		"""
		# find segment in which integration should occur
		for i, seg in enumerate(self.model[integration.chr]):
			if seg['origin'] != 'host':
				continue
			if integration.pos >= seg['coords'][0] and integration.pos <= seg['coords'][1]:
				break
		
		t0 = time.time()
		# remove this element from the list
		seg = self.model[integration.chr].pop(i)
		host_start = seg['coords'][0]
		host_stop = seg['coords'][1]
		left_overlap_bases = integration.junc_props['n_junc'][0] if integration.junc_props['junc_types'][0] == 'overlap' else 0
		right_overlap_bases = integration.junc_props['n_junc'][1] if integration.junc_props['junc_types'][1] == 'overlap' else 0
		overlap_bases = left_overlap_bases + right_overlap_bases
		
		# create host segment before this integration and add to list
		host = {'origin' : 'host', 'seq_name' : integration.chr, 'ori' : '+'}
		
		# if the integration had a left overlap, we need to trim these bases from the host
		# note that int.pos is always to the left of any overlaps
		# so we don't need to consider them here
		host['coords'] = (host_start, integration.pos)
		
		assert host['coords'][1] > host['coords'][0]
		self.model[integration.chr].insert(i, host)
		i += 1

		# if we have ambiguous bases at the left junction, add these to the list too
		assert len(integration.junc_props['junc_bases'][0]) == integration.junc_props['n_junc'][0]
		if integration.junc_props['junc_types'][0] in ['gap', 'overlap']:
			# features common to both ambiguous types
			ambig = {'origin': 'ambig',
					 'ori' : "+"}
			ambig['bases'] = integration.junc_props['junc_bases'][0]
			
			# gap features
			if integration.junc_props['junc_types'][0] == 'gap':
				ambig['seq_name'] = 'gap'
				ambig['coords'] = (0, integration.junc_props['n_junc'][0])
				
			# overlap features
			elif integration.junc_props['junc_types'][0] == 'overlap':
				ambig['seq_name'] = integration.chr
				ambig['coords'] = (integration.pos, integration.pos + left_overlap_bases)
				assert str(self.host[integration.chr][ambig['coords'][0]:ambig['coords'][1]].seq).lower() == integration.junc_props['junc_bases'][0].lower()
			else:
				raise ValueError(f"unrecgonised integration type: {integration.junc_props[0]}")
			assert ambig['coords'][1] > ambig['coords'][0]
			self.model[integration.chr].insert(i, ambig)
			i += 1

		# add each piece of the viral chunk too
		for j in range(len(integration.chunk.pieces)):
			virus = {'origin': 'virus', 
					 'coords': (integration.chunk.pieces[j][0], integration.chunk.pieces[j][1]),
					 'ori' : integration.chunk.oris[j],
					 'seq_name' : integration.chunk.virus}
			assert virus['coords'][1] > virus['coords'][0]
			self.model[integration.chr].insert(i, virus)
			i += 1
			
		t0 = time.time()	
		# if we have ambiguous bases at the right junction, add these
		assert len(integration.junc_props['junc_bases'][1]) == integration.junc_props['n_junc'][1]
		if integration.junc_props['junc_types'][1] in ['gap', 'overlap']:
			ambig = {'origin': 'ambig',
					 'bases' : integration.junc_props['junc_bases'][1],
					 'ori' : "+"}
					 
			if integration.junc_props['junc_types'][1] == 'gap':
				ambig['seq_name'] = 'gap'
				ambig['coords'] = (0, integration.junc_props['n_junc'][1])
				
			# overlap features
			elif integration.junc_props['junc_types'][1] == 'overlap':
				ambig['seq_name'] = integration.chr
				
				# if the left junction was also an overlap, we need to account for this in the coordinates
				if integration.junc_props['junc_types'][0] == 'overlap':
					start = integration.pos + left_overlap_bases
					stop = start + right_overlap_bases
				else:
					start = integration.pos
					stop = start + right_overlap_bases
				
				ambig['coords'] = (start, stop)
				assert str(self.host[integration.chr][ambig['coords'][0]:ambig['coords'][1]].seq).lower() == integration.junc_props['junc_bases'][1].lower()
			else:
				raise ValueError(f"unrecgonised integration type: {integration.junc_props[0]}")
			assert ambig['coords'][1] > ambig['coords'][0]
			self.model[integration.chr].insert(i, ambig)
			i += 1
		
			
		# finally, add second portion of host
		host = {'origin': 'host', 'seq_name': integration.chr, 'ori': '+'}
		
		
		# accounting for bases deleted from the host and overlaps
		# note that integration.pos is always to the left of any overlapped bases, so here is where we need to account for them
		host['coords'] =  (integration.pos + integration.junc_props['host_del_len'] + overlap_bases, host_stop)

		assert host['coords'][1] > host['coords'][0]
		
		self.model[integration.chr].insert(i, host)
		i += 1		
		
		if self.model_check is True:
			self.__check_model()
		
	
	def __check_model(self):
		"""
		check model is valid by checking various properties
		"""
		n_ints = 0
		next_int = True
		t0 = time.time()
		for chr in self.model.keys():
			host_pos = 0
			for seg in self.model[chr]:
				assert seg['coords'][1] > seg['coords'][0]
				assert seg['origin'] in ('host', 'virus', 'ambig')
				assert 'seq_name' in seg
				if seg['origin'] == 'host':
					next_int = True
					assert seg['seq_name'] == chr
					assert seg['ori'] == '+'
					# check that host position is only increasing
					assert seg['coords'][0] >= host_pos
					host_pos = seg['coords'][1]
				elif seg['origin'] == 'virus':
					if next_int is True:
						n_ints += 1
						next_int = False	
					assert seg['seq_name'] in list(self.virus.keys())
				elif seg['origin'] == 'ambig':
					assert 'bases' in seg
					if seg['seq_name'] != 'gap':
						assert seg['seq_name'] in list(self.host.keys())
						host_bases = str(self.host[chr][seg['coords'][0]:seg['coords'][1]].seq).lower()
						seg_bases = seg['bases'].lower()
						assert host_bases == seg_bases
		assert n_ints == len(self)
		t1 = time.time()
		print(f"checked model validity in {t1-t0}s")	
	
	def __save_fasta(self, filename, append = False):
		"""
		Save host fasta with integrated viral bases to a fasta file
		"""
		assert isinstance(append, bool)
		if append is True:
			write_type = "a"
		if append is False:
			write_type = "w+"
			
		with open(filename, write_type) as handle:
			# loop over chromosomes
			for chr in self.host.keys():
			
				# print chromosome name
				handle.write(f">{chr}\n")
								
				# loop over entries in the model for this chromosome
				for entry in self.model[chr]:
					start = entry['coords'][0]
					stop = entry['coords'][1]
					
					# if host
					if entry['origin'] == 'host':
						handle.write(str(self.host[chr].seq[start:stop]))
				
					# if ambiguous write bases - note that overlapped bases have been trimmed from host and virus
					# so we're ok to write them here
					elif entry['origin'] == 'ambig':
						handle.write(entry['bases'])
							
					# if viral
					elif entry['origin'] == 'virus':
						virus = entry['seq_name']
						if entry['ori'] == '+':
							handle.write(str(self.virus[virus].seq[start:stop]))
						elif entry['ori'] == '-':
							handle.write(str(self.virus[virus].seq[start:stop].reverse_complement()))
						else:
							raise ValueError(f"unregconised orientation {entry['ori']} in {entry}")
					
					else:
						raise ValueError(f"unrecgonised model feature on chr {chr}: {entry}")
				

				handle.write("\n")
		
	def __save_info(self, filename):
		"""
		Output the following info for each integration (one integration per line) into a tab-separated file:
		Note that all co-orindates are 0-based
		 - chr: host chromosome on which integration occurs
		 - hPos: 0-based position in the original host fasta (before adding variation with simuG) at which integration occurs 
		 - hPos_input_fasta: 0-based position in the input host fasta at which integration occurs
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
		#pp.pprint(self.model)
		
		# dictionary will keep track of the number of bases previously integrated on each chromosome
		previous_ints = {key:0 for key in self.host.keys()}
		deleted_bases = {key:0 for key in self.host.keys()}
		
		self.__check_model()
		self.sort()
			
		with open(filename, "w", newline='') as csvfile:
			intwriter = csv.writer(csvfile, delimiter = '\t')
			intwriter.writerow(['id', 'chr', 'hPos', 'hPos_input_fasta', 'leftStart', 'leftStop',  'rightStart', 'rightStop', 'hDeleted', 'virus', 'vBreakpoints', 'vOris', 'juncTypes', 'juncBases', 'juncLengths', 'whole', 'rearrangement', 'deletion', 'n_swaps', 'n_delete', 'in_indel'])
			for i, integration in enumerate(self):
				assert integration.pos is not None

				# to calculate hPos, we need to account for variants introduced by simuG
				h_pos = integration.pos
				
				# find all relevant indels (on the same chromosome and before this integration position)
				if self.indel is not None:
					rel_indels = [row for row in self.indel if integration.chr == row['ref_chr']]
					rel_indels = [row for row in rel_indels if integration.pos >= int(row['sim_end'])]
				
					# is this integration in an indel?
					in_indel = any([integration.pos >= int(row['sim_start']) and integration.pos <= int(row['sim_end']) for row in rel_indels])
					
					# net number of bases added/deleted
					net_len = sum([len(row['sim_allele']) - len(row['ref_allele']) for row in rel_indels])
					
					# adjust position
					h_pos += net_len
				
				# calculate start and stop position for this integration			
				left_start = integration.pos + previous_ints[integration.chr] - deleted_bases[integration.chr]
				left_stop = left_start + integration.junc_props['n_junc'][0] 
				
				right_start = left_stop + len(integration.chunk.bases)
				right_stop = right_start + integration.junc_props['n_junc'][1] 
				
				# update previous_ints - total integrated bases
				# are the integrated viral bases, and the bases in the gaps
				previous_ints[integration.chr] += len(integration.chunk.bases)
				if integration.junc_props['junc_types'][0] == 'gap':
					previous_ints[integration.chr] += integration.junc_props['n_junc'][0]
				if integration.junc_props['junc_types'][1] == 'gap':
					previous_ints[integration.chr] += integration.junc_props['n_junc'][1]
				
				# update deleted_bases
				deleted_bases[integration.chr] += integration.junc_props['host_del_len']

				# format lists into comma-separated strings
				breakpoints = ";".join([str(i) for i in integration.chunk.pieces])
				oris = ','.join(integration.chunk.oris)
				junc_types = ",".join(integration.junc_props['junc_types'])
				junc_bases = ",".join(integration.junc_props['junc_bases'])
				junc_lengths = ",".join([str(i) for i in integration.junc_props['n_junc']])
				
				# write a row for this integration
				intwriter.writerow([integration.id,
									integration.chr, 
									h_pos, 
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
									integration.chunk.chunk_props['is_whole'],
									integration.chunk.chunk_props['n_swaps'] > 0,
									integration.chunk.chunk_props['n_delete'] > 0,
									integration.chunk.chunk_props['n_swaps'],
									integration.chunk.chunk_props['n_delete'],
									in_indel])

	def __str__(self):
		return f"Viral integrations object with {len(self)} integrations of viral sequences {list(self.virus.keys())} into host chromosomes {list(self.host.keys())}"

	def __repr__(self):
		return f"Object of type Integrations with properties {self}"	
	
class Episomes(Integrations):
	"""
	Episomes may be added to the output fasta to mimic contamination of sample with purely viral sequences
	this class stores a list of ViralChunk objects which make up the contaminating viral sequences
	
	This class is intended to be used by the Integrations class
	
	Since Integrations and Episomes use some similar methods, this class inherits from Integrations
	in order to avoid duplication
	"""
	def __init__(self, virus, rng, probs, max_attempts = 50, min_len = None):
		"""
		initialises an empty Episomes list, storing the properties common to all episomes
		"""
		# checks for inputs
		assert isinstance(virus, dict)
		assert isinstance(probs, dict)
		assert isinstance(max_attempts, int)
		assert isinstance(min_len, int) or min_len is None
		if min_len is None:
			min_len = 1
		else:
			assert min_len > 0
		
	
		# assign properties common to all episomes
		self.virus = virus
		self.probs = probs
		self.rng = rng
		self.max_attempts = max_attempts
		self.min_len = min_len
		
		# default values for probs
		self.set_probs_default('p_whole', default_p_whole)
		self.set_probs_default('p_rearrange', default_p_rearrange)
		self.set_probs_default('p_delete', default_p_delete)
		self.set_probs_default('lambda_split', default_lambda_split)
		
		# checks on assigned variables
		assert 'p_whole' in probs
		assert 'p_rearrange' in probs
		assert 'p_delete' in probs	
		assert 'lambda_split' in probs
			
		# check that probabilities are between 0 and 1
		self.check_prob(self.probs['p_whole'], 'p_whole')
		#self.check_prob(self.probs['p_rearrange'] + self.probs['p_delete'], 'the sum of p_rearrange and p_delete')	
		self.check_prob(self.probs['p_rearrange'], 'p_rearrange')
		self.check_prob(self.probs['p_delete'], 'p_delete')
		self.check_prob(self.probs['p_overlap'] + self.probs['p_gap'], 'the sum of p_overlap and p_gap')

		# when adding episomes, we will want to keep track of the number of episomes for each virus
		self.virus_counts = {virus: 0 for virus in self.virus.keys()}
	
	def add_episome(self):
		"""
		get a viral chunk and add to self
		"""
		
		# call functions that randomly set properties of this integrations
		chunk_props = self.set_chunk_properties()
		assert len(chunk_props) != 2
		
		# get a viral chunk
		chunk = ViralChunk(self.virus, self.rng, chunk_props)
		
		# check for valid chunk
		if chunk.pieces is not None:
			# add id to chunk
			chunk.id = f"{chunk.virus}__{self.virus_counts[chunk.virus]}"
			self.virus_counts[chunk.virus] += 1
			# append to self
			self.append(chunk)
			return True
		
		return False
	
	def __save_fasta(self, filename, append = True):
		"""
		save each ViralChunk as a separate sequence in an output fasta
		"""
		
		assert isinstance(append, bool)
		if append is True:
			write_type = "a"
		if append is False:
			write_type = "w+"
		
		# virus counts will keep track of the number of episomes
		# for each virus
		virus_counts = {virus: 0 for virus in self.virus.keys()}
		
		# open file for writing
		with open(filename, write_type) as handle:
		
			for chunk in self:
				
				# write name of virus and current count as an id
				handle.write(f">{chunk.id}\n")
				
				virus_counts[chunk.virus] += 1
				
				# write viral bases
				handle.write(f"{chunk.bases}\n")
		
	def __save_info(self, filename):
		"""
		save chunk.chunk_properties into tab-separated format
		fields to write:
		 - id: a unique identifier for each episome
		 - virus: virus from which this episome comes from
		 - start: coordinate (in virus) of left-most base in chunk
		 - stop: coordinate (in virus) of right-most base in chunk
		 - pieces: breakpoints of each piece of this chunk
		 - oris: orientation of each piece of this chunk
		 - is_whole: if this episome was originally the whole virus
		 - is_rearrange: if this episome was rearranged
		 - is_deletion: if this episome contains a deletion
		 - n_swaps: number of swaps made when rearranging
		 - n_delete: number of pieces deleted from chunk
		
		"""
		assert len(self) >= 0
		
		# define header
		
		header = ["id", "virus", "start", "stop", "pieces", "oris", "is_whole", "is_rearrange", "is_deletion", "n_swaps", "n_delete"]
		
		with open(filename, "w", newline='') as csvfile:
			epiwriter = csv.writer(csvfile, delimiter = '\t')
			epiwriter.writerow(header)
			
			for chunk in self:
				epiwriter.writerow([chunk.id,
								 chunk.virus,
								 chunk.start,
								 chunk.stop,
								 chunk.pieces,
								 chunk.oris,
								 chunk.chunk_props['is_whole'],
								 chunk.chunk_props['n_swaps'] > 0,
								 chunk.chunk_props['n_swaps'],
								 chunk.chunk_props['n_delete'] > 0,
								 chunk.chunk_props['n_delete']])
			
	def __str__(self):
		return f"Episomes object with {len(self)} episomes drawn from viral sequences {list(self.virus.keys())}"

	def __repr__(self):
		return f"Object of type Episomes with properties {self}"
		
class Integration(dict):
	"""
	Class to store the properties of an individual integration.  If properly instantiated, it stores
	the properties of the junctions either side of the integrated bases, and a ViralChunk object (self.chunk)
	which stores the properties of the integrated bases themselves.
	
	This class is intended to be used by the Integrations class, which is essentially a list of Integration objects
	"""
	def __init__(self, host, virus, model, rng, int_id, chunk_props, junc_props, min_sep):
		"""
		Function to initialise Integration object
		portionType is 'whole' or 'portion' - the part of the virus that has been inserted
		overlapType is two-member tuple of 'clean', 'gap' or 'overlap' - defines the junction at each end of the integration
		
		when initialising an Integration, need to provide:
		 - host (as a dict of SeqRecord objects or similar)
		 - virus (as a dict of SeqRecord ojbects or similar)
		 - model - self.model of Integrations object - specifies where existing integrations have occured
		 - rng (a numpy.random.Generator object for setting random properties)
		 - chunk_props (a dict of properties for initialising ViralChunk object)
		 - junc_props (a dict of properties defining the junctions of this integration)
		 	- junc_types = iterable of length 2 specifying type of left and right junctions (one of gap, overlap, clean)
			- n_junc = iterable of length 2 specifying length of left and right junction
		
		objects of this class have the following properties:
		self.id - an integer unique to this integrations
		self.chr - string: host chromosome on which integration should be located
		self.chunk - ViralChunk object which is integrated
		self.chunk_props - properties of the integration chunk that were specified as input. dict with fields:
			- is_whole: boolean specifying if the ViralChunk is whole (if false, chunk is just a portion)
			- n_fragments: number of fragments into which ViralChunk should be split
			- n_delete: number of fragments to delete from ViralChunk (should always leave at least two fragments after deletion)
			- n_swaps: number of swaps to make when rearranging ViralChunk
		self.junc_props - properties of the integration chunk that were specified as input
			- junc_types = iterable of length 2 specifying type of left and right junctions (one of gap, overlap, clean)
			- n_junc = iterable of length 2 specifying length of left and right junction
			- junc_bases = iterable of length 2 specifying bases involved in each junction
			- host_del = boolean specifying if there should be a deletion from the host
			- host_del_len = number of bases to be deleted from the host at this integration site
			
		if anything went wrong with initialisation, self.pos is set to None - check this to make sure integration
		is valid
		
		after assigning a 
		"""
		# assign chunk_props and junc_props to self
		self.chunk_props = chunk_props
		self.junc_props = junc_props

		# check inputs	
		assert isinstance(virus, dict)
		
		assert isinstance(model, dict)
		for key in host.keys():
			assert key in model
		
		assert isinstance(chunk_props, dict) # leave most of the checking to ViralChunk
		
		assert isinstance(junc_props, dict)
		assert len(junc_props) == 4
		
		assert 'junc_types' in junc_props
		assert len(self.junc_props['junc_types']) == 2
		assert all([i in ['clean', 'gap', 'overlap'] for i in self.junc_props['junc_types']])
		
		assert 'n_junc' in junc_props
		assert len(self.junc_props['n_junc']) == 2
		assert all([isinstance(i, int) for i in self.junc_props['n_junc']])
		assert all([i >=0 for i in self.junc_props['n_junc']])
		
		assert 'host_del' in junc_props
		assert isinstance(junc_props['host_del'], bool)
		
		assert 'host_del_len' in junc_props
		assert isinstance(junc_props['host_del_len'], int)
		assert junc_props['host_del_len'] >= 0
		if junc_props['host_del'] is True:
			assert junc_props['host_del_len'] > 0
		if junc_props['host_del'] is False:
			assert junc_props['host_del_len'] == 0
		
		assert search_length_overlap > 0 and isinstance(search_length_overlap, int)
		
		assert isinstance(min_sep, int)
		assert min_sep > 1
			
		# set parameters that won't be changed
		self.search_length_overlap = search_length_overlap
		self.id = int_id
	
		# get random chromosome on which to do integration
		self.chr = str(rng.choice(list(host.keys())))
		
		# set minimum length for viral chunk - longer than the number of bases involved in the junction
		if self.chunk_props['min_len'] is None:
			self.chunk_props['min_len'] = self.n_overlap_bases() + 1
		if self.chunk_props['min_len'] < self.n_overlap_bases() + 1:
			self.chunk_props['min_len'] = self.n_overlap_bases() + 1
		
		# get viral chunk
		self.chunk = ViralChunk(virus, 
								rng, 
								self.chunk_props)
		# if specified properties (chunk_props) are incompatible with the initialised chunk
		# self.chunk.pieces will be None
		if self.chunk.pieces is None:
			self.pos = None
			return
		
		
		# set the number of bases in overlaps
		self.junc_props['n_junc'] = [junc_props['n_junc'][0], junc_props['n_junc'][1]]
		# but overwrite in the case of a clean junction
		if self.junc_props['junc_types'][0] == 'clean':
			self.junc_props['n_junc'][0] = 0
		if self.junc_props['junc_types'][1] == 'clean':
			self.junc_props['n_junc'][1] = 0

		
		# number of bases in overlaps must be less than the length of the integrated chunk
		if self.n_overlap_bases() >= len(self.chunk.bases):
			self.pos = None
			return
		
		# set bases belonging to junction
		self.junc_props['junc_bases'] = (self.get_junc_bases(rng, 'left'), self.get_junc_bases(rng, 'right'))

		
		# get a position at which to integrate
		pos_success = self.get_int_position(host[self.chr].seq, rng, model, min_sep)
		if pos_success is False:
			self.pos = None
			return

		
		# number of bases deleted from host chromosome
		if junc_props['host_del'] is False:
			self.host_deleted = 0 
		else:
			self.host_deleted = junc_props['host_del_len']
			
			# but only delete up to the length of the segment in which this integration occurs, 
			# so that we don't delete any integrations as well
			for seg in model[self.chr]:
		
				# skip viral and ambiguous segments
				if seg['seq_name'] != self.chr:
					continue
					
				# find if this is the segment in which the integration occurs
				if not self.has_overlap(self.pos, self.pos, seg['coords'][0], seg['coords'][1]):
					continue
			
				# are we trying to delete past the end of the segment?		
				if self.pos + self.host_deleted + self.n_overlap_bases() >= (seg['coords'][1] - min_sep):
					self.host_deleted = seg['coords'][1] - self.pos - min_sep - self.n_overlap_bases()
					self.junc_props['host_del_len'] = self.host_deleted
					if self.host_deleted < 0:
						self.pos = None
						return
				break
			
		# double check for valid chunk
		assert self.chunk.bases == self.chunk.get_bases(virus)
		assert 'junc_bases' in self.junc_props
		assert len(self.junc_props['junc_bases']) == 2
		assert len(self.junc_props['junc_bases'][0]) == self.junc_props['n_junc'][0]
		assert len(self.junc_props['junc_bases'][1]) == self.junc_props['n_junc'][1]
		assert all([len(i) == 2 for i in self.chunk.pieces])
		assert len(self.chunk.pieces) == len(self.chunk.oris)
		assert all([piece[1] > piece[0] for piece in self.chunk.pieces])

		
	def n_overlap_bases(self):
		"""
		Get the total number of bases in overlaps
		"""
		assert len(self.junc_props['n_junc']) == 2
		assert len(self.junc_props['junc_types']) == 2	
		n = 0
		if self.junc_props['junc_types'][0] == 'overlap':
			n += self.junc_props['n_junc'][0]
		if self.junc_props['junc_types'][1] == 'overlap':
			n += self.junc_props['n_junc'][1]
		return n
		
	def get_random_position(self, model, rng, min_sep):
		"""
		based on current model, get a random position that is available for integration
		(that is, doesn't already have an integration)
		
		do this by getting all the host parts of the model, and choosing a random one (weighted by their length)
		and then choosing a random position from within this part
		"""

		# get a list of host coordinates
		host_coords = [segment['coords'] for segment in model[self.chr] if segment['origin'] == 'host']
		
		# TODO - enforce minimum separation
		host_coords = [(coords[0] + min_sep, coords[1] - min_sep) for coords in host_coords]
		
		# ensure lengths are positive
		host_coords = [coord for coord in host_coords if (coord[1] - coord[0]) > 0]
		
		# get lengths of each range to weight choosing a part
		lengths = [coord[1] - coord[0] for coord in host_coords]
		lengths = [length/sum(lengths) for length in lengths]
		
		# get a random part, weighted by lengths
		part = rng.choice(host_coords, p = lengths)
		
		# get a random postion from this part
		return int(rng.choice(range(part[0], part[1])))	
		
	def get_int_position(self, chr, rng, model, min_sep):
		"""
		get a position at which to perform the integration
		the position depends on the overlap type at each side of the overlap
		if both sides are 'clean' or 'gap', position can be random
		'overlap' junctions place constraints on the integration location because the need
		to be placed where there are overlaps
		"""

		# if at both ends the junction is either 'clean' or 'gap', just get a random positon
		if all([True if (i in ['clean', 'gap']) else False for i in self.junc_props['junc_types']]):
			self.pos = self.get_random_position(model, rng, min_sep)
			return True
			
		# if overlap at both ends, look for both overlaps in host chromosome next to each other
		elif self.junc_props['junc_types'][0] == 'overlap' and self.junc_props['junc_types'][1] == 'overlap':
		
			# make string with both overlap bases 
			self.pos = self.find_overlap_bases(self.junc_props['junc_bases'][0] + self.junc_props['junc_bases'][1], chr, rng, model, min_sep)
			
			# check for unsuccessful find
			if self.pos == -1:
				self.pos = None
				return False
			
			## need to remove overlapped bases from viral chunk, and adjust chunk start, breakpoints and oris
			self.delete_left_bases(self.junc_props['n_junc'][0])
			self.delete_right_bases(self.junc_props['n_junc'][1])
			
			return True
		
		# if one end is an overlap, find those bases in the host chromosome
		# left overlap
		elif self.junc_props['junc_types'][0] == 'overlap':
			
			# find position with overlap at which to do overlap
			self.pos = self.find_overlap_bases(self.junc_props['junc_bases'][0], chr, rng, model, min_sep)
			# check for unsuccessful find
			if self.pos == -1:
				self.pos = None
				return False
			
			## need to remove overlapped bases from viral chunk, and adjust chunk start, breakpoints and oris
			self.delete_left_bases(self.junc_props['n_junc'][0])
			
			return True

		# right overlap
		elif self.junc_props['junc_types'][1] == 'overlap':

			# find position with overlap at which to do overlap
			self.pos = self.find_overlap_bases(self.junc_props['junc_bases'][1], chr, rng, model, min_sep)
			# check for unsuccessful find
			if self.pos == -1:
				self.pos = None
				return False
			
			## need to remove overlapped bases from viral chunk, and adjust chunk start, breakpoints and oris
			self.delete_right_bases(self.junc_props['n_junc'][1])
			
			return True
			
		else:
			raise ValueError(f"junction types {self.junc_props['junc_types']} are not implemented yet")	
			
		return False		
	
	def find_overlap_bases(self, overlap_bases, chr, rng, model, min_sep):
		"""
		find bases from an overlap in the host chromosome
		"""
		# get position around which to search
		start = self.get_random_position(model, rng, min_sep)
		
		# get start and stop positions for bases in host chromosome to search for overlaps
		stop = start + self.search_length_overlap
		
		# make sure that we aren't searching after the end of the chromosome
		if stop > len(chr):
			stop = len(chr)
			
		# find overlapping bases in host segments in model
		for seg in model[self.chr]:
			# check that this segment comes from the host
			if seg['origin'] != 'host':
				continue
				
			# check that this segment overlaps with desired start and stop from above			
			if not self.has_overlap(start, stop, seg['coords'][0], seg['coords'][1]):
				continue
				
			# look for desired overlap bases in this segment
			search_start = seg['coords'][0] + min_sep
			search_stop = seg['coords'][1] - min_sep
			
			if search_stop <= search_start:
				continue
			
			if seg['ori'] == "+": 
				found = re.search(overlap_bases, str(chr[search_start:search_stop]), re.IGNORECASE)
				if found:
					return found.span()[0] + search_start
			else:
				# not implemented - we assume no flipping of host chromosome segments
				raise ValueError("host chromosome segment {seg} has a negative orientation!")
				
		
		# check for unsuccessful find
		return -1
		
	def has_overlap(self, start_1, stop_1, start_2, stop_2):
		"""
		check to see if two intervals, specified by start_1 and stop_1; and start_2 and stop_2, overlap
		"""
		# check inputs
		assert isinstance(start_1, int)
		assert isinstance(start_2, int)
		assert isinstance(stop_1, int)
		assert isinstance(stop_2, int)
		assert start_1 >= 0
		assert start_2 >= 0
		assert stop_1 >= 0
		assert stop_2 >= 0
		assert start_1 <= stop_1
		assert start_2 <= stop_2
		
		# interval 1 start and stop to the left of interval 2
		if (start_1 < start_2) and (stop_1 < start_2):
			return False
			
		# interval 2 start and stop are to the right of interval 2
		if (start_1 > stop_2) and (stop_1 > stop_2):
			return False
		
		# otherwise they must overlap
		return True
			
	def delete_left_bases(self, n):
		"""
		delete bases on the left after adding an overlap - need to adjust chunk bases, oris and breakpoints and start
		"""
		
		assert self.junc_props['junc_types'][0] == 'overlap'
		
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
			# if we're left with a piece of length 0, flag this piece for deletion
			if self.chunk.pieces[i][0] == self.chunk.pieces[i][1]:
				to_delete.append(i)
				i += 1
				
		# remove chunks that we want to delete
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
		
		assert self.junc_props['junc_types'][1] == 'overlap'
		
		# check we're not trying to delete more bases than there are in the chunk
		assert n < len(self.chunk.bases)
		
		# adjust stop
		self.chunk.stop -= n
		
		# adjust bases
		self.chunk.bases = self.chunk.bases[:-n]
		
		# adjust breakpoints 		
		deleted_bases = 0
		to_delete = []
		i = len(self.chunk.pieces) - 1 # start at the last piece in the chunk
		while deleted_bases < n:
			assert i >= 0
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
			# if we're left with a piece of length 0, flag this piece for deletion			
			if self.chunk.pieces[i][0] == self.chunk.pieces[i][1]:
				to_delete.append(i)
				i -= 1
				 
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
		assert len(self.junc_props['junc_types']) == 2
		assert len(self.junc_props['n_junc']) == 2
		assert self.junc_props['junc_types'][0] in ['gap', 'overlap', 'clean']
		assert self.junc_props['junc_types'][1] in ['gap', 'overlap', 'clean']

		
		if side == 'left':
			n_bases = self.junc_props['n_junc'][0]
			# no bases in a clean junction
			if self.junc_props['junc_types'][0] == 'clean':
				return ""
			# random bases in a gap
			elif self.junc_props['junc_types'][0] == 'gap':
				return self.get_n_random_bases(rng, n_bases)
			# first n bases of viral chunk in an overlap
			elif self.junc_props['junc_types'][0] == 'overlap':
				return self.chunk.bases[:n_bases]
			else:
				raise ValueError(f"unrecgonised type: {self.junc_props['junc_types'][0]}")
		elif side == 'right':
			n_bases = self.junc_props['n_junc'][1]
			if self.junc_props['junc_types'][1] == "clean":
				return ""
			elif self.junc_props['junc_types'][1] == 'gap':
				return self.get_n_random_bases(rng, n_bases)
			elif self.junc_props['junc_types'][1] == 'overlap':
				return self.chunk.bases[-n_bases:]
			else:
				raise ValueError(f"unrecgonised type: {self.junc_props['junc_types'][1]}")
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
	Class to store properties of an integrated chunk of virus.  
	Intended to be used by the Integrations and Episomes classes
	"""

	def __init__(self, virus,  rng, chunk_props):
		"""
		function to get a chunk of a virus
		virus is the dictionary of seqrecords imported using biopython
		
		
		desired properties of chunk are input as a dict with the following attributes:
			- is_whole: boolean specifying if the ViralChunk is whole (if false, chunk is just a portion)
			- n_fragments: number of fragments into which ViralChunk should be split
			- n_delete: number of fragments to delete from ViralChunk (should always leave at least two fragments after deletion)
			- n_swaps: number of swaps to make when rearranging ViralChunk
			- min_len: the minimum length of this chunk (integer greater than 1)
		
		the bases attribute of a ViralChunk consist of only the bases that are unique to the virus. 
		So in the case of an Integration of a ViralChunk with a 'overlap' type junction,
		the bases, breakpoints and oris attributes are re-assigned to remove the overlapped bases
		
		"""
		# check inputs
		assert isinstance(virus, dict)
		assert len(virus) > 0

		assert isinstance(chunk_props, dict)
		assert len(chunk_props) == 5
		
		assert 'is_whole' in chunk_props
		assert isinstance(chunk_props['is_whole'], bool)
		
		assert 'n_fragments' in chunk_props
		assert isinstance(chunk_props['n_fragments'], int)
		assert chunk_props['n_fragments'] > 0
		
		assert 'n_delete' in chunk_props
		assert isinstance(chunk_props['n_delete'], int)
		assert chunk_props['n_delete'] >= 0
		
		assert 'n_swaps' in chunk_props
		assert isinstance(chunk_props['n_delete'], int)
		assert chunk_props['n_delete'] >= 0
		
		assert 'min_len' in chunk_props
		assert isinstance(chunk_props['min_len'], int) or chunk_props['min_len'] is None
		if chunk_props['min_len'] is None:
			chunk_props['min_len'] = 1
		assert chunk_props['min_len'] > 0
		
		# check that the minimum length specified is longer than all the viruses
		# otherwise we might fail
		if not all([chunk_props['min_len'] <= len(vir.seq) for vir in virus.values()]):
			raise ValueError("minimum length must be longer than all the Viruses")
		
		# get a random virus to integrate
		self.virus = str(rng.choice(list(virus.keys())))
		self.chunk_props = chunk_props
		
		# if the length of the virus is equal to min_len, is_whole must be true
		if len(virus[self.virus]) == chunk_props['min_len']:
			chunk_props['is_whole'] = True
		
		if self.chunk_props['is_whole'] is True:
			self.start = 0
			self.stop = len(virus[self.virus])
		elif self.chunk_props['is_whole'] is False:
			self.start = int(rng.integers(low = 0, high = len(virus[self.virus].seq) - chunk_props['min_len']))
			self.stop = int(rng.integers(low = self.start + chunk_props['min_len'], high = len(virus[self.virus].seq)))
		else:
			raise ValueError("self.chunk_props['is_whole'] must be either True or False")
			
		# breakpoints are the start and stop coordinates of pieces of the virus that have been 
		# integrated
		
		# set breakpoints
		self.pieces = [[self.start, self.stop]]
		
		self.oris = str(rng.choice(('+', '-')))
		
		# do a deletion if applicable
		
		if self.chunk_props['n_delete'] > 0:
			self.__delete(rng)
			# if something went wrong, breakpoints will be None
			if self.pieces is None:
				return
				
		if self.chunk_props['n_swaps'] > 0:
			self.__rearrange(rng)
			# if something went wrong, breakpoints will be None
			if self.pieces is None:
				return
		
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
		
	def __split_into_pieces(self, rng):
		"""
		get random, unique breakpoints to divide a viral chunk into pieces
		there must be at least min_breakpoints, which results in min_breakpoints + 1 pieces 
		this is a list of coordinates, not tuples (unlike self.pieces)
		"""
		# shouldn't do this if this chunk has already been split
		assert len(self.pieces) == 1
		assert len(self.oris) == 1
		
		# check that we're not trying to divide a chunk into more pieces than there are bases
		if self.chunk_props['n_fragments'] >= self.stop - self.start:
			self.pieces = None
			return
		
		# get the number of pieces to divide into
		num_breakpoints = self.chunk_props['n_fragments'] - 1
				
		# get random breakpoints from within this chunk
		breakpoints = rng.choice(range(self.start + 1, self.stop - 1), size = num_breakpoints, replace = False)
				
		# set self.pieces
		breakpoints = [self.start] + sorted(breakpoints) + [self.stop]
		self.pieces = [[int(breakpoints[i]), int(breakpoints[i+1])] for i in range(len(breakpoints) - 1)]
		
		# set self.oris
		self.oris = [self.oris[0]] * len(self.pieces)
		return	
		
	def __swap_orientations(self, breakpoint, side):
		"""
		Given a breakpoint, swap all of the orientations (+ to - or vice versa) for all of the pieces
		on the left or right of this breakpoint
		"""
		if side == 'left':
			for i, ori in enumerate(self.oris[:breakpoint]):
				if ori == "+":
					self.oris[i] = "-"
				else:
					self.oris[i] = "+"
		else:
			for i, ori in enumerate(self.oris[breakpoint:]):
				if ori == "+":
					self.oris[i] = "-"
				else:
					self.oris[i] = "+"		
	
	def __delete(self, rng):
		"""
		Divide a viral chunk up into multiple pieces
		and remove one of those pieces
		"""
		
		# deletions are always performed first, so the chunk should not have been split yet
		assert len(self.pieces) == 1
		
		# want to have at least two pieces left
		assert self.chunk_props['n_fragments'] - self.chunk_props['n_delete'] >= 2 

		# split chunk into at n_fragments pieces
		self.__split_into_pieces(rng)
		if self.pieces is None:
			return
		assert len(self.pieces) == self.chunk_props['n_fragments']
		
		# decide which portions to delete
		i_delete = rng.choice(range(1, len(self.pieces) - 1), self.chunk_props['n_delete'], replace=False)
		
		# do deletion
		self.pieces = [piece for i, piece in enumerate(self.pieces) if i not in i_delete]
		self.oris = [ori for i, ori in enumerate(self.oris) if i not in i_delete]	
		
		assert len(self.pieces) == self.chunk_props['n_fragments'] - self.chunk_props['n_delete']
		
	def __rearrange(self, rng):
		"""
		Divide a viral chunk up into multiple pieces
		and randomise their order and orientiations
		"""
		# split the chunk if it hasn't already been split
		if len(self.pieces) == 1:
			# split chunk into at least three pieces
			self.__split_into_pieces(rng)
			if self.pieces is None:
				return
			assert len(self.pieces) == self.chunk_props['n_fragments']
		else:
			assert len(self.pieces) > 1
		
		# if we only have two pieces, we should only do one swap
		# so that we don't end up back with the same fragment
		# there are other ways to end up with the same fragment after swaps
		# but don't worry about them for now - TODO
		if len(self.pieces) == 2:
			self.chunk_props['n_swaps'] = 1
		
		for i in range(self.chunk_props['n_swaps']):
			# pick a point about which to swap
			if 1 == len(self.pieces) - 1:
				i_swap = 1
			else:
				i_swap = rng.choice(range(1, len(self.pieces) - 1))
				
			# swap everything to the left of this position with everything on the right
			self.pieces = self.pieces[i_swap:] + self.pieces[:i_swap]
			
			# 50 % chance of swapping the orientations of all the pieces for each side
			if bool(rng.choice((True, False))) is True:
				self.__swap_orientations(i_swap, 'left')
			if bool(rng.choice((True, False))) is True:
				self.__swap_orientations(i_swap, 'right')

	def __str__(self):
		return f"Viral chunk of virus {self.virus} ({self.start}, {self.stop}) and orientations {self.oris}"
		
	def __repr__(self):
		return f"Object of type ViralChunk with properties {self}"

if __name__ == "__main__":
	main(argv[1:])
