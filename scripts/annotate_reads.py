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
from os import path, SEEK_CUR
from scipy.stats import norm
from collections import defaultdict
import argparse
import csv
import pysam

import pdb

genome_length = int(3e9)

def main(argv):
	#get arguments
	parser = argparse.ArgumentParser(description='simulate viral integrations')
	parser.add_argument('--sim-info', help='information from simulation', required = True, type=str)
	parser.add_argument('--sim-sam', help='sam file from read simulation', required = True, type=str)
	parser.add_argument('--soft-threshold', help='threshold for number of bases that must be mapped when finding integrations in soft-clipped reads', type=int, default=20)
	parser.add_argument('--discordant-threshold', help='threshold for number of bases that may be unmapped when finding integrations in discordant read pairs', type=int, default=20)
	parser.add_argument('--mean-frag-len', help='mean framgement length used when simulating reads', type=int, default=500)
	parser.add_argument('--sd-frag-len', help='standard devation of framgement length used when simulating reads', type=int, default=30)
	parser.add_argument('--window-frac', help='fraction of the distribution of fragment sizes to use when searching for discordant read pairs', type=float, default=0.99)
	parser.add_argument('--output', help='output file', required=False, default='annotated.tsv')
	args = parser.parse_args()
	
	# read in bam/sam file
	samfile = pysam.AlignmentFile(args.sim_sam)
	
	# iterate over integrations in info file and pull out reads crossing each one
	# the sam/bam file has the same coordinates as the info file
	# so just use the coordinates of the left and right junctions from this files
	with open(args.sim_info, newline='') as info_file, open(args.output, 'w', newline='') as output:
		
		# use csv to read and write files
		reader = csv.DictReader(info_file, delimiter = '\t')
		
		writer_fieldnames = list(reader.fieldnames) + ['left_chimeric', 'right_chimeric', 'left_discord', 'right_discord', 'multiple_discord', 'fake_discord', 'discordant_and_chimeric']
		
		writer = csv.DictWriter(output, delimiter = '\t', fieldnames = writer_fieldnames)
		writer.writeheader()
		
		
		# set window size for looking for discordant pairs and 
		# looking for multiple ints with the same discordant pair
		window_width = window_size(args.mean_frag_len, args.sd_frag_len, args.window_frac)
		
		# create a buffer of integrations that all fall within this window width 
		buffer = [next(reader)]
		while True:

			buffer, next_buffer = get_ints_in_window(buffer, window_width, reader, info_file)
			
			# find reads crossing each integration in the buffer
			buffer = find_reads_crossing_ints(buffer, samfile, args, window_width)
		
			# check for read pairs that cross multiple junctions
			buffer = find_multiple_discordant(buffer)
			
			# check for read pairs that are discordant at one integration and chimeric at another integration
			# assume that integrations are independent (even though they aren't in simulation), so these are not detectable
			buffer = find_multiple_chimeric_discordant(buffer)
			
			# write rows in buffer
			for row in buffer:
				writer.writerow(row)
			
			buffer = next_buffer
			
			# check if we've reached the end of the file
			if len(buffer) == 0:
				break
		
	samfile.close()
	
def find_multiple_chimeric_discordant(buffer):
	"""
	Assume that in real data, integrations are independent, so two integrations do not
	occur close enough to each other in the same cell that we can detect them.  This may
	or may not be the case, but the number of these kinds of pairs of integrations is assumed
	to be small enough that they're not worth worrying about. 
	
	If a read pair has one read that is chimeric, but it's also discordant about the same
	or a different integration, only the chimeric read will be detected 
	(and it won't be flagged as a discordant read pair). Assuming that
	these events are a tiny minority in real data, here don't flag the discordant pairs
	as something that should be detected (but instead add them to the 'discordant_and_chimeric'
	category)
	"""	
	
	to_delete = {}
	for i, row in enumerate(buffer):
	
		# for each row, find reads that are both chimeric
		# and discordant for any other integration
		left_discord = row['left_discord'].split(';')
		right_discord =  row['right_discord'].split(';')
		to_delete[i] = {'left' : [], 'right' : []}
		
		# check for these reads in other rows
		for row_2 in buffer:
			# skip if this is the same as the one we're looking at
			if row_2 == row:
				continue
				
			# get read names for chimeric reads for this row
			chimeric = row_2['left_chimeric'].split(";")
			chimeric += row_2['right_chimeric'].split(";")
			chimeric = [read[:-2] for read in chimeric]
			
			# check if chimeric reads in row_2 also occur in the list
			# of discordant reads in row
			for read in chimeric:
				if read in left_discord:
					to_delete[i]['left'].append(read)
				if read in right_discord:
					to_delete[i]['right'].append(read)
			
	# remove reads in to_delete
	for row_num in to_delete.keys():
		
		# remove reads from left_discord
		left_discord = buffer[row_num]['left_discord'].split(';')
		left_discord = [read for read in left_discord if read not in to_delete[row_num]['left']]
		buffer[row_num]['left_discord'] = ";".join(left_discord)
		
		# remove reads from right_discord
		right_discord = buffer[row_num]['right_discord'].split(';')
		right_discord = [read for read in right_discord if read not in to_delete[row_num]['right']]
		buffer[row_num]['right_discord'] = ";".join(right_discord)		
	
		# removed reads go to 'discordant_and_chimeric'
		buffer[row_num]['discordant_and_chimeric'] = [f"{read}_left" for read in to_delete[row_num]['left']]
		buffer[row_num]['discordant_and_chimeric'] += [f"{read}_right" for read in to_delete[row_num]['right']]
		buffer[row_num]['discordant_and_chimeric'] = ";".join(buffer[row_num]['discordant_and_chimeric'])
	
	return buffer

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
			if read == '':
				continue
			if read in seen.keys():
				# if this is the first double, add the entry from seen
				if read not in multiples:
					multiples[read] = []
					multiples[read].append(seen[read])
				
				# add entry to multiples for this time we found the read
				multiples[read].append({'int_id' : row['id'], 'side' : 'left' , 'buffer_row' : i})
				
			else:
				seen[read] = {'int_id' : row['id'], 'side' : 'left', 'buffer_row' : i}
		# check right side
		for read in row['right_discord'].split(';'):
			if read == '':
				continue
			if read in seen.keys():
				# if this is the first double, add the entry from seen
				if read not in multiples:
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
		return len(multiples[read])
	else:
		return 0
	
def find_reads_crossing_ints(buffer, samfile, args, window_width):
	"""
	find reads crossing the integration site, and add them to the row
	reads crossing the integration site can be chimeric about the left or right junction
	or discordant about the left or right junction
	if they're discordant about both junctions, they actually go from host sequence to host sequence
	and	therefore aren't actually discordant
	"""
	for row in buffer:
		
		# get information about integration site location
		chr = row['chr']
		left_start = int(row['leftStart']) 
		left_stop = int(row['leftStop'])
		right_start = int(row['rightStart']) 
		right_stop = int(row['rightStop'])
	
		assert left_start >= 0 and right_start >= 0

		# find chimeric reads
		left_chimeric = get_chimeric(chr, left_start, left_stop, samfile, args.soft_threshold, buffer)
		right_chimeric = get_chimeric(chr, right_start, right_stop, samfile, args.soft_threshold, buffer)
	
		row['left_chimeric'] = ";".join(left_chimeric)
		row['right_chimeric'] = ";".join(right_chimeric)
	
		# find discordant read pairs	
		left_discord = get_discordant(chr, left_start, left_stop, samfile, args.discordant_threshold, window_width, buffer)
		right_discord = get_discordant(chr, right_start, right_stop, samfile, args.discordant_threshold, window_width, buffer)
	
		# if a read is both chimeric and discordant, chimeric takes priority 
		# (but it shouldn't be both if the clipping threshold is the same for both read types)
		left_chimeric, left_discord = remove_chimeric_from_discord(left_chimeric, left_discord)
		right_chimeric, right_discord = remove_chimeric_from_discord(right_chimeric, right_discord)
			
		row['left_discord'] = ";".join(left_discord)
		row['right_discord'] = ";".join(right_discord)
	
	return buffer	
			
def get_ints_in_window(buffer, window_width, reader, reader_handle):
	"""
	get a list of integrations (dicts corresponding to one row from the int-info file)
	that are within window_width of each other
	"""	
		
	assert len(buffer) == 1
	
	# get position and chromosome from buffer
	start = int(buffer[0]['rightStop'])
	chr = buffer[0]['chr']
	
	# get more integraions until we're not in the window anymore or we're at the end of the file
	# to avoid having to go back a line, save the first line not added to 'next_buffer'
	while True:	
		
		# get next row
		try:
			row = next(reader)
		except StopIteration:
			next_buffer = []
			break
		
		# compare previous integration with this integration
		prev_start = start
		prev_chr = chr
		start =  int(row['leftStart'])
		chr = row['chr']
		
		# check if next row is a window width away
		if (start - prev_start < window_width) and prev_chr == chr:
			# add to buffer
			buffer.append(row)
			start = int(row['rightStop'])
		else:
			# don't add the row but this will be the start of the next buffer
			next_buffer = [row]
			break

	return buffer, next_buffer	

def get_discordant(chr, start, stop, samfile, threshold, window_width, buffer):
	"""
	Get any discordant read pairs which cross an integration site
	In other words, get pairs where one mate is mapped on the host side, and the other on the virus side
	A read is considered mapped if it has at most threshold (default 20) unmapped bases
	
	This includes any pairs where the integration site falls within threshold (default 20) bases 
	of the end of the read (for a read mapped on the left of the integration), 
	or within threshold bases of the start of the read for a read on the mapped on the right
	of the integration
	
	Avoid an exhaustive search by extracting only the reads in a window around the integration site
	Set this window based on the mean length and standard deviation of fragment size used in simulation
	and the fraction of the fragment length distribution we want to cover.  Find discordant
	pairs by finding pairs for which an integration (start, stop) falls between the 20th
	base from the end of the left read and the 20th base of the right read
	
	Current criteria in discordant.pl is that a pair must have one read mapped and the other
	unmapped to be considered a discordant read pair.  To be considered 'mapped', a read must
	have 20 or fewer unmapped bases.  To be considered 'unmapped', a read must have
	less than 20 unmapped bases
	
	Note that a read pair might have one integration fall between read1 and read2, but read1 or read2
	might also cross a second integration.  This pair is therefore both discordant about one integration, and
	also one member of the pair is chimeric  A consequence of this is that one read maps only to vector or host, 
	but the other maps to both.  This discordant pair cannot currently be detected, since
	the pipeline currently detects discordant read-pairs only if one read is mapped to vector and
	not host, and vice versa for the other read.  However, this pair really is evidence
	for integration at the first site (as a discordant read-pair), so include it in the
	output as such.
	"""
	
	reads = []
	
	# extract read pairs in the desired window
	# pysam numbering is 0-based, with the only exception being the region string in the fetch() and pileup() methods. 
	window_start = start - window_width
	if window_start < 1:
		
		window_start = 1
	
	window_stop = stop + window_width
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
# 			
# 			# TODO - need to decide if integrations should be excluded on the basis
# 			# below (if it's not the case that one read is mapped to host and the other to virus)
# 			# these events can't be currently detected in the pipeline, but are real evidence
# 			# of integration, so for now include them in the output
# 
# 			r1_mapped = get_mapped_ref(read1, buffer, threshold)
# 			r2_mapped = get_mapped_ref(read2, buffer, threshold)
# 			
# 			assert r1_mapped['host'] or r1_mapped['virus']
# 			assert r2_mapped['host'] or r2_mapped['virus']
# 			
# 			if r1_mapped['host'] != r2_mapped['host'] and r1_mapped['virus'] != r2_mapped['virus']:
# 				reads.append(read1.qname)
	
	return reads
	
def get_mapped_ref(read, buffer, threshold):
	"""
	figure out if each read in this read pair will be mapped to host or vector/virus
	returns a dict with the keys 'host' and 'virus',
	and values True or False depending on if the read is mapped to either or not
	"""
	
	assert read.is_unmapped is False
	
	read_mapped = {'host':False, 'virus':False, 'int' : []}
	
	# get first and last position to which read is mapped in reference
	first = read.get_reference_positions()[0]
	last =  read.get_reference_positions()[-1]
	
	prev_start = 0
	for row in buffer:
		# figure out if we need to include the ambiguous bases or not
		if row['juncTypes'].split(',')[0] == 'gap':
			left_host_junc = int(row['leftStart'])
			left_virus_junc = int(row['leftStop'])
		else:
			left_host_junc = int(row['leftStop'])
			left_virus_junc = int(row['leftStart'])
		if row['juncTypes'].split(',')[1] == 'gap':
			right_virus_junc = int(row['rightStart'])
			right_host_junc = int(row['rightStop'])
		else:
			right_virus_junc = int(row['rightStop'])
			right_host_junc = int(row['rightStart'])	

		# if read is between the start of the chromosome and the start of the integration
		if intersects(first, last, prev_start, left_host_junc):
			# check we have at least threshold bases mapped
			if check_threshold(read, prev_start, left_host_junc, threshold):
				read_mapped['host']  = True
		# if read is between the left and right junctions
		if intersects(first, last, left_virus_junc, right_virus_junc):
			if check_threshold(read, left_virus_junc, right_virus_junc, threshold):
				read_mapped['virus'] = True
				read_mapped['int'].append(row['id'])
		prev_start = right_host_junc
		
		# if the prev_start is to the right of the position of the last mapped base in the read
		# then we don't need to continue
		if prev_start > last:
			break
		
				
	# check if read is between the end of the last integration and the end of the chromosome
	# don't know what the end of the chromosome is, so just use a large number (genome_length)
	if intersects(first, last, prev_start, genome_length):
		if check_threshold(read, prev_start, genome_length, threshold):
			read_mapped['host'] = True
	
	return read_mapped	
	
def get_chimeric(chr, start, stop, samfile, threshold, buffer):
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
		mapped = get_mapped_ref(read, buffer, threshold)
		assert (mapped['host'] or mapped['virus'])
		
		if mapped['host'] is True and mapped['virus'] is True:
		
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
	check that there are least threshold bases that map to an interval (defined by start and stop)
	"""
	rstart = read.get_reference_positions()[0]
	# need to account for 0-based numbering (stop already has this accounted for)
	rstop = read.get_reference_positions()[-1] + 1
	assert intersects(rstart, rstop, start, stop)
	
	if (rstop - start) < threshold:
		return False
	if (stop - rstart) < threshold:
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
	get first position in reference that is at least threshold bases from the 
	start or end of the read and isn't None
	
	if side is 'left', return the 0-based position after threshold mapped bases from the 
	start of the read
	if side is 'right', return the 0-based position before threshold mapped bases from the
	end of the read
	"""
	assert isinstance(threshold, int)
	assert threshold >= 0
	assert threshold <= read.qlen
	assert side in ['left', 'right']
		
	aligned_pairs = read.get_aligned_pairs()
	
	if side == "left":
		# if we want to look on the right hand side, we need to look backwards 
		# through the aligned pairs
		aligned_pairs = list(reversed(aligned_pairs))
		# python numbering is zero-based and half-open
		threshold -= 1 

	
	# iterate over bases in aligned_pairs, starting from the threshold value
	for pair in aligned_pairs[threshold:]:
		# check the base is aligned
		if (pair[1] is not None):
			return pair[1]
		
def within(start1, stop1, start2, stop2):
	"""
	compare two intervals, each with a start and stop value
	return true if the first interval is completely within in the second
	
	use half-open intervals, so [8, 8) is not within [5, 8) 
	within(8, 8, 5, 8) => False
	within(6, 6, 5, 8) => True
	within(5, 8, 5, 8) => True
	within(4, 6, 5, 8) => False
	within(5, 9, 5, 8) => False
	"""	
	assert start1 <= stop1
	assert start2 <= stop2
	assert isinstance(start1, int)
	assert isinstance(start2, int)
	assert isinstance(stop1, int)
	assert isinstance(stop2, int)
	
	if start1 >= start2 and start1 < stop2:
		if stop1 - 1 >= start2 and stop1 - 1 < stop2:
			return True
	return False
		
def intersects(start1, stop1, start2, stop2):
	"""
	compare two intervals, each with a start and stop value
	return true if the first interval is intersects the second
	
	use half-open intervals, so [2, 3) and [3, 4) don't intersect
	"""	
	assert start1 <= stop1
	assert start2 <= stop2
	assert isinstance(start1, int)
	assert isinstance(start2, int)
	assert isinstance(stop1, int)
	assert isinstance(stop2, int)
	
	# if first interval is completely to the left of the second interval
	if stop1 <= start2:
		return False
	# if first interval is completely to the right of the second interval
	if start1 >= stop2:
		return False
		
	# otherwise, they intersect
	return True	

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
	to_delete = []
	chimeric = [read[:-2] for read in chimeric]
	for read in discord:
		if read in chimeric:
			print(f"  removed a read that was in both chimeric and discord: {read}")
			to_delete.append(read)
	
	# remove reads flagged for deletion
	discord = [read for read in discord if read not in to_delete]
	
	return chimeric, discord

if __name__ == "__main__":
	main(argv[1:])
