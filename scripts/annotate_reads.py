# annotate info file with information about which reads belong to which integration

from sys import argv
from os import path
from scipy.stats import norm
from collections import defaultdict
import argparse
import csv
import pysam


import pdb
import re

def main(argv):
	#get arguments
	parser = argparse.ArgumentParser(description='simulate viral integrations')
	parser.add_argument('--sim-info', help='information from simulation', required = True, type=str)
	parser.add_argument('--sim-sam', help='sam file from read simulation', required = True, type=str)
	parser.add_argument('--soft-threshold', help='threshold for amount of read that must be mapped when finding integrations in soft-clipped reads', type=int, default=20)
	parser.add_argument('--mean-frag-len', help='mean framgement length used when simulating reads', type=int, default=500)
	parser.add_argument('--sd-frag-len', help='standard devation of framgement length used when simulating reads', type=int, default=30)
	parser.add_argument('--window-frac', help='fraction of the distribution of fragment sizes to use when searching for discordant read pairs', type=float, default=0.99)
	parser.add_argument('--output', help='output file', required=False, default='results.tsv')
	args = parser.parse_args()
	
	# read in bam/sam file
	samfile = pysam.AlignmentFile(args.sim_sam)
	
	# iterate over integrations in info file and pull out reads crossing each one
	# the sam/bam file has the same coordinates as the info file
	# so just use the coordinates of the left and right junctions from this file
	with open(args.sim_info, newline='') as info_file, open(args.output, 'w', newline='') as output:
		
		# use csv to read and write files
		reader = csv.DictReader(info_file, delimiter = '\t')
		
		writer_fieldnames = list(reader.fieldnames) + ['left_chimeric', 'right_chimeric', 'left_discord', 'right_discord', 'both_discord']
		
		writer = csv.DictWriter(output, delimiter = '\t', fieldnames = writer_fieldnames)
		writer.writeheader()
		
		for row in reader:
			# pull out reads that cross the left junction
			chr = row['chr']
			left_start = int(row['leftStart']) 
			left_stop = int(row['leftStop'])
			right_start = int(row['rightStart']) 
			right_stop = int(row['rightStop'])
			
			assert left_start >= 0 and right_start >= 0
			
			print(f"finding reads for integration {row['id']}")

			# find chimeric reads
			left_chimeric = get_chimeric(chr, left_start, left_stop, samfile, args.soft_threshold)
			right_chimeric = get_chimeric(chr, right_start, right_stop, samfile, args.soft_threshold)
			
			row['left_chimeric'] = ";".join(left_chimeric)
			row['right_chimeric'] = ";".join(right_chimeric)
			
			# find discordant read pairs
			window_width = window_size(args.mean_frag_len, args.sd_frag_len, args.window_frac)
			left_discord = get_discordant(chr, left_start, left_stop, samfile, args.soft_threshold, window_width)
			right_discord = get_discordant(chr, right_start, right_stop, samfile, args.soft_threshold, window_width)
			
			# if a read is both chimeric and discordant, chimeric takes priority

			left_chimeric, left_discord = remove_chimeric_from_discord(left_chimeric, left_discord)
			right_chimeric, right_discord = remove_chimeric_from_discord(right_chimeric, right_discord)
					
			# if a read crosses both the left and right boundaries, it will just be human-human and won't be picked up
			both_discord = []
			for read in left_discord:
				if read in right_discord:
					both_discord.append(read)
					left_discord.remove(read)
					right_discord.remove(read)
			
			row['left_discord'] = ";".join(left_discord)
			row['right_discord'] = ";".join(right_discord)
			row['both_discord'] = ";".join(both_discord)
			

			writer.writerow(row)
		
	samfile.close()

def get_discordant(chr, start, stop, samfile, threshold, window_width):
	"""
	Get any discordant read pairs which cross an integration site
	In other words, get pairs where one mate is mapped on the host side, and the other on the virus side
	
	This includes any pairs where the integration site falls within threshold bases of the end
	of the read, on the side of the read that is closest to the mate
	
	Avoid an exhaustive search by extracting only the reads in a window around the integration site
	Set this window based on the mean length and standard deviation of fragment size used in simulation
	and the fraction of the fragment length distribution we want to cover.
	"""
	
	reads = []
	
	# extract read pairs in the desired window
	# pysam numbering is 0-based, with the only exception being the region string in the fetch() and pileup() methods. 
	window_start = start - int(round(window_width / 2)) - 1
	if window_start < 1:
		
		window_start = 1
	
	window_stop = stop + int(round(window_width / 2))
	if window_stop > samfile.get_reference_length(chr):
		window_stop = samfile.get_reference_length(chr) 
	if window_stop == window_start:
		window_stop += 1
	
	
	for read1, read2 in read_pair_generator(samfile, f"{chr}:{window_start}-{window_stop}"):
		# check mate is mapped
		if read1.is_unmapped or read2.is_unmapped:
			continue
		# check reference for this read is the same as mate
		if read1.reference_name != read2.reference_name:
			continue
		# if this read is forward, mate must be reverse and vice versa
		if (read1.is_reverse == read2.is_reverse):
			continue
		
		# if the integration site falls between left_boundary and right_boundary
		# (which are coordinates within the reference)
		# this pair crosses the integration site
		if read1.is_reverse is False:
			left_boundary = get_boundary(read1, threshold, side = "left")
			right_boundary =  get_boundary(read2, threshold, side = "right")
		else:
			left_boundary = get_boundary(read2, threshold, side = "left")
			right_boundary =  get_boundary(read1, threshold, side = "right")

		assert left_boundary is not None
		assert right_boundary is not None
		assert left_boundary < right_boundary
		
		if within(start, stop, left_boundary, right_boundary):
			reads.append(read1.qname)
	
	return reads
	
def get_chimeric(chr, start, stop, samfile, threshold):
	"""
	find reads that cross an interval defined as chr:start-stop in samfile
	the interval must be at least threshold bases from the start and end of the read
	"""
	reads = []
	
	# get reads that cross interval
	# pysam numbering is 0-based, with the only exception being the region string in the fetch() and pileup() methods. 
	# The same is true for any coordinates passed to the samtools command utilities directly, such as pysam.fetch().
	for read in samfile.fetch(chr, start, stop + 1):
		# check that interval is at least threshold bases from either end of the read
		if check_threshold(read, start, stop, threshold) is False:
			continue
			
			
		reads.append(read.query_name + read_num(read))
		
	return reads
	
def read_num(read):
	"""
	return '/1' if read is R1, '/2' if read is R2, or empty string otherwise
	"""
	if read.is_read1 is True:
		return "/1"
	elif read.is_read2 is True:
		return "/2"
	else:
		return ""
		
def check_threshold(read, start, stop, threshold):
	""""
	check that there are least threshold bases between an integration site (defined by start and stop)
	and the start and end of the read
	"""
	
	if (start - read.get_reference_positions()[0]) < threshold:
		return False
	if (read.get_reference_positions()[-1] - stop) < threshold:
		return False
		
	return True
	
def window_size(mean_frag_len, sd_frag_len, window_frac):
	"""
	to avoid exhaustive search for discordant reads, 
	work out window size for extracting reads when finding discordant read pairs
	based on mean and standard deviation of fragment length distribution,
	and fraction of this distribution we want to cover
	
	For example, if the mean fragment length is 500 bp, the standard deviation is 30 bp, and we want 
	to cover 0.99 of this distribution, the window size would be 570 bp, which accounts for 99% of 
	the fragment sizes in this distribution (one-tailed)
	
	"""
	
	# get one-tailed value which contains window_frac of fragments
	upper = norm.ppf(window_frac, loc = mean_frag_len, scale = sd_frag_len)
	
	return int(round(upper))
	
def get_boundary(read, threshold, side = "left"):
	"""
	get first position in reference that is at least threshold bases from the end of the read
	and isn't None
	"""
	assert isinstance(threshold, int)
	assert threshold >= 0
	assert side in ['left', 'right']
	if side == "left":
		aligned_pairs = reversed(read.get_aligned_pairs())
	else:
		aligned_pairs = read.get_aligned_pairs()
	
	for pair in aligned_pairs:
		# for left side, we want to check if we're at least threshold bases
		# before the end of the read
		if side == "left":
			if ((read.qlen - pair[0] - 1) >= threshold) and (pair[1] is not None):
				return pair[1]
		# for the right side, we want to be at least threshold bases i
		else:
			if (pair[0] >= threshold) and (pair[1] is not None):
				return pair[1]
		
def within(start1, stop1, start2, stop2):
	"""
	compare two intervals, each with a start and stop value
	return true if the first interval is completely encompassed by the second
	"""	
	assert start1 <= stop1
	assert start2 <= stop2
	
	if start1 >= start2 and stop1 <= stop2:
		return True
	else:
		return False

def read_pair_generator(bam, region_string=None):
    """
    https://www.biostars.org/p/306041/
    Generate read pairs in a BAM file or within a region string.
    Reads are added to read_dict until a pair is found.
    """
    read_dict = defaultdict(lambda: [None, None])
    for read in bam.fetch(region=region_string):
        if read.is_secondary or read.is_supplementary:
            continue
        qname = read.query_name
        if qname not in read_dict:
            if read.is_read1:
                read_dict[qname][0] = read
            else:
                read_dict[qname][1] = read
        else:
            if read.is_read1:
                yield read, read_dict[qname][1]
            else:
                yield read_dict[qname][0], read
            del read_dict[qname]
            
def remove_chimeric_from_discord(chimeric, discord):
	"""
	check for read ids that are in both chimeric and discord - remove them from discord if they're in both
	"""
	# chimeric reads have /1 or /2 added
	chimeric = [read[:-2] for read in chimeric]
	for i, read in enumerate(discord):
		if read in chimeric:
			print(f"  removed a read that was in both chimeric and discord: {discord[i]}")
			del discord[i]
			
	return chimeric, discord

if __name__ == "__main__":
	main(argv[1:])