# annotate info file with information about which reads belong to which integration

# for each integration, assign reads that are informative about that integration to one of 
# the following categories:
# - left_chimeric: reads that are chimeric at the left boundary of the integration
# - right_chimeric: reads that are chimeric at the right boundary of the integration
# - left_discord: read pairs where one read maps in the host before the integration, and one read maps in the viral sequence
# - right_discord: read pairs where one read maps in the viral sequence, and one read maps in the host sequence after the integration
# - multiple_discord: read pairs where one read maps to virus and one to host, but span more than one integration 

# there is also the following category of reads that are indirect evidence of integration
# - fake_discord: read pairs where both reads map to either host or vector, but span more than one integration 
#				  in other words, read pairs where one both reads map to either host or vector, but are interrupted by a sequence of the other type
#				  for example, a pair where one read maps to the virus of one integration, and the other to the virus of the next integration


from sys import argv
from os import path
from scipy.stats import norm
from collections import defaultdict
import argparse
import csv
import pysam

import pdb

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
	# so just use the coordinates of the left and right junctions from this files
	with open(args.sim_info, newline='') as info_file, open(args.output, 'w', newline='') as output:
		
		# use csv to read and write files
		reader = csv.DictReader(info_file, delimiter = '\t')
		
		writer_fieldnames = list(reader.fieldnames) + ['left_chimeric', 'right_chimeric', 'left_discord', 'right_discord', 'multiple_discord', 'fake_discord']
		
		writer = csv.DictWriter(output, delimiter = '\t', fieldnames = writer_fieldnames)
		writer.writeheader()
		
		
		# set window size for looking for discodant pairs and 
		# looking for multiple ints with the same discordant pair
		window_width = window_size(args.mean_frag_len, args.sd_frag_len, args.window_frac)
		
		# create a buffer of integrations that all fall within this window width 
		while True:

			buffer = get_ints_in_window(window_width, reader)
			
			# check if we've reached the end of the file
			if len(buffer) == 0:
				break
		
			# find reads crossing each integration in the buffer
			buffer = [find_reads_crossing_ints(row, samfile, args, window_width)
						for row in buffer]
		
			# check for read pairs that cross multiple junctions
			buffer = find_multiple_discordant(buffer)
			
			# write rows in buffer
			for row in buffer:
				writer.writerow(row)
		
	samfile.close()
	
def find_multiple_discordant(buffer):
	"""
	a read might be discordant about two or more integration sites
	for example, if a read pair crosses the right border of one integration and the left
	border of a close by integration, both read pairs will be mapped to the virus
	and therefore this pair won't be detected as an integration
	
	deal with discordant pairs that cross multiple integrations in one of two ways:
		1. if a pair is discordant about the right side of one integration and the left side
		of a close by integration, then both reads are mapped to the virus and this read should
		be removed from both.  Keep track of these in a different category ('fake_discord')
		2. if a pair crosses right side of one integration, and both sides of a nearby integration,
		then one read maps to virus (in the first integration), and the other maps to host (after the second)
		and therefore this pair should be indicated to be an integration.  However, if we want to score based
		on integration position, we need to decide which integration this read 'belongs' to.  Keep
		these reads in a seperate category ('multiple_discord') for all the integrations that they cross,
		to reflect the fact that the 'belong' to multiple integrations
	"""
	
	# keep track of which reads we've already seen
	seen = dict()
	
	# for reads that we find at multiple junctions, 
	# keep track of which junctions so we can work out what to do with that read
	multiples = dict()
	
	for i, row in enumerate(buffer):
		# check reads found around left side
		for read in row['left_discord'].split(';'):
			if read in seen.keys():
				# find out how many times we've already seen this read ID
				n = times_already_found(read, multiples)
				
				# if this is the first double, add the entry from seen
				if n == 0:
					multiples[read] = []
					multiples[read].append(seen[read])
				
				# add entry to multiples for this time we found the read
				multiples[read].append({'int_id' : row['id'], 'side' : 'left' , 'buffer_row' : i})
				
			else:
				seen[read] = {'int_id' : row['id'], 'side' : 'left', 'buffer_row' : i}
		# check right side
		for read in row['right_discord'].split(';'):
			if read in seen.keys():
				# find out how many times we've already seen this read ID
				n = times_already_found(read, multiples)
				
				# if this is the first double, add the entry from seen
				if n == 0:
					multiples[read] = []
					multiples[read].append(seen[read])
				
				# add entry to multiples for this time we found the read
				multiples[read].append({'int_id' : row['id'], 'side' : 'right', 'buffer_row' : i})
			else:
				seen[read] = {'int_id' : row['id'], 'side' : 'right', 'buffer_row' : i}
		
		# add extra keys to dicts 
		row['multiple_discord'] = []
		row['fake_discord'] = []
			
	# deal with reads crossing multiple integrations
	for read in multiples:
	
		# if the side of the first junction matches the side of the last junction
		# this is a 'multiple_discord'
		sides = [find['side'] for find in multiples[read]]
		if sides[0] == sides[-1]:
			new_type = 'multiple_discord'
		else:
			new_type = 'fake_discord'
	
		for find in multiples[read]:
			
			# find row and side
			buffer_row = find['buffer_row']
			old_type = f"{find['side']}_discord"
			
			# remove from list of reads of old type
			reads = buffer[buffer_row][old_type].split(";")
			reads.remove(read)
			buffer[buffer_row][old_type] = ";".join(reads)
			
			# add to list of new type
			buffer[buffer_row][new_type].append(f"{read}_{find['side']}")

	# for each row, join lists of 'multiple_disord' and 'fake_discord'
	for row in buffer:
		row['multiple_discord'] = ";".join(row['multiple_discord'])
		row['fake_discord'] = ";".join(row['fake_discord'])
	
	return buffer

def times_already_found(read, multiples):
	"""
	check how many times a read has already been found when checking for discordant about multiple integration sides
	"""	
	if read in multiples:
		return max(multiples[read].keys())
	else:
		return 0
	
def find_reads_crossing_ints(row, samfile, args, window_width):
	"""
	find reads crossing the integration site, and add them to the row
	reads crossing the integration site can be chimeric about the left or right junction
	or discordant about the left or right junction
	if they're discordant about both junctions, they actually go from host sequence to host sequence
	and	therefore aren't actually discordant
	"""
	
	# get information about integration site location
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
	
	left_discord = get_discordant(chr, left_start, left_stop, samfile, args.soft_threshold, window_width)
	right_discord = get_discordant(chr, right_start, right_stop, samfile, args.soft_threshold, window_width)
	
	# if a read is both chimeric and discordant, chimeric takes priority 
	# (but it shouldn't be both if the clipping threshold is the same for both read types)
	left_chimeric, left_discord = remove_chimeric_from_discord(left_chimeric, left_discord)
	right_chimeric, right_discord = remove_chimeric_from_discord(right_chimeric, right_discord)
			
	row['left_discord'] = ";".join(left_discord)
	row['right_discord'] = ";".join(right_discord)
	
	return row	
			
def get_ints_in_window(window_width, reader):
	"""
	get a list of integrations (dicts corresponding to one row from the int-info file)
	that are within window_width of each other
	"""	
	
	buffer = []
		
	# get first integration, unless we're at the end of the file
	try:
		first_row = next(reader)
	except StopIteration:
		return buffer
	
	first_stop = int(first_row['rightStop'])
	last_start = int(first_row['rightStop'])
	first_chr = first_row['chr']
	last_chr = first_row['chr']
	buffer.append(first_row)
	
	# get more integraions until we're not in the window anymore or we're at the end of the file
	while (last_start - first_stop < window_width) and (first_chr == last_chr):
	
		try:
			row = next(reader)
			last_start = int(row['leftStart'])
			last_chr = row['chr']
			buffer.append(row)
		except StopIteration:
			break	
		
		
	return buffer		

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
			left_boundary = get_boundary(read1, threshold, side = "right")
			right_boundary =  get_boundary(read2, threshold, side = "left")
		else:
			left_boundary = get_boundary(read2, threshold, side = "right")
			right_boundary =  get_boundary(read1, threshold, side = "left")
		
		# if left_boundary is greater than right_boundary, reads overlap
		if left_boundary >= right_boundary:
			continue
		
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
	assert threshold <= read.qlen
	assert side in ['left', 'right']
		
	aligned_pairs = read.get_aligned_pairs()
	
	# if we want to look on the right hand side, we need to look backwards through the aligned pairs
	if side == "right":
		aligned_pairs = list(reversed(aligned_pairs))
	
	# iterate over bases in aligned_pairs, starting from the threshold value
	for pair in aligned_pairs[threshold:]:
		# check the base is aligned
		if (pair[1] is not None):
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
