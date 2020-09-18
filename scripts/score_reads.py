#!/usr/bin/env python3

# Score simulated vs detected integrations on a per-read basis

# score each read as a true positive, true negative, false positive or false negative
# also for reads which are found in simulation and analysis results, check if properties are the same

# simplest scoring is that a read is:
#### true positive if it crosses at least one integration (ie, is annotated in one line in the simulation information file)
####	and is also found in at least one line of the analysis file

#### true negative if it doesn't cross any integrations (ie, isn't annotated in any lines in the simulation information file)
####	and is also not found in any lines of the analysis file

#### false negative if it doesn't cross any integrations (ie, is annotated in one line in the simulation information file)
####	but is not found in any lines of the analysis file

#### false positive if it doesn't cross any integrations (ie, isn't annotated in any lines in the simulation information file)
####	but is found in at least one line of the analysis file

# however, might also wish to consider not only if read is found in the simulation and analysis results
# but also if it is found for the correct integration in both results

# to speed things up, process a buffer of reads in parallel and use callback function to collate
# and write results to file

from sys import argv
from os import cpu_count
import argparse
import csv
import pdb
import pprint
import multiprocessing as mp
import functools
import time


n_print = 10000

def main(argv):
	#get arguments
	parser = argparse.ArgumentParser(description='simulate viral integrations')
	parser.add_argument('--sim-info', help='information from simulation with reads annotated', required = True, type=str)
	parser.add_argument('--analysis-info', help='information from analysis of simulated reads', required = True, type=str)
	parser.add_argument('--sim-sam', help='sam file from ART with simulated reads (must not be binary!)', required=True, type=str)
	parser.add_argument('--wiggle-room', help='allow this number of bases wiggle room when deciding if sim and analysis reads overlap, and start and stop are in the same place', type=int, default=5)
	parser.add_argument('--output', help='output file for results', required=False, default='results.tsv')
	parser.add_argument('--output-summary', help='file for one-line output summary', required=False, default='results-summary.tsv')
	parser.add_argument('--threads', help='number of threads to use', required=False, type=int)
	parser.add_argument('--polyidus', help='analysis of simulated reads was performed with polyidus', action='store_true')
	parser.add_argument('--vifi', help='analysis of simulated reads was performed with ViFi', action='store_true')
	
	args = parser.parse_args()
	
	# check that we haven't specified polyidus and vifi
	if args.vifi and args.polyidus:
		raise ValueError("Data can only come from either ViFi or Polyidus (not both)")
	
	# read simulated and pipeline information into memory
	print(f"opening simulated information: {args.sim_info}")
	sim_info = read_csv(args.sim_info)

	print(f"opening analysis information: {args.analysis_info}")
	if args.polyidus is True:
		analysis_info = read_polyidus_csv(args.analysis_info)
	elif args.vifi is True:
		analysis_info = read_vifi_file(args.analysis_info)
	else:
		analysis_info = read_csv(args.analysis_info)
	
	print(f"opening sam file: {args.sim_sam}")
	with open(args.sim_sam) as samfile:

		# keep count of true and false positives and negatives
		read_scores = {'tp':0, 'tn':0, 'fp':0, 'fn':0}
		score_types = ['found_score', 'host_score', 'virus_score']
		read_scores = {score_type : {'chimeric' : dict(read_scores), 'discord' : dict(read_scores)} for score_type in score_types}
	
		# create DictWriter object for output
		output_header =  ['readID', 'intID'] + score_types + ['found', 'n_found',
					'side', 'correct_side', 'type', 'correct_type', 'correct_host_chr',
					'host_start_dist', 'host_stop_dist', 'host_coords_overlap', 'correct_virus', 'virus_coords_overlap', 'virus_start_dist',
					'virus_stop_dist', 'ambig_diff']
	
		# create queues
		line_queue = mp.Queue()
		result_queue = mp.Queue()
 		
 		# create pool of workers
		if args.threads is None:
			args.threads = cpu_count()
		workers = args.threads - 1
		if workers < 1:
			workers = 1
 			
 		# create output process
		print(f"using {workers} workers")
		output_p = mp.Process(target = write_results, args = (output_header, result_queue, read_scores, workers, args))
		output_p.start()
 		
		workers = [mp.Process(target = process_lines, 
 					args =  (line_queue, result_queue, sim_info, analysis_info, args))
 					for i in range(workers)]
		[worker.start() for worker in workers]
 		
 		# iterate over lines in samfile		
		[line_queue.put(line) for line in samfile]
 		
 		# let the worker processes know that there's no more lines
		[line_queue.put(None) for i in range(args.threads)]		
 		
		[worker.join() for worker in workers]	
 		
		output_p.join()

		# for processing serially
#		with open(args.output, "w", newline = '') as outfile:
#			output_writer = csv.DictWriter(outfile, delimiter = '\t', 
#										fieldnames = output_header)
#			output_writer.writeheader()
#			read_count = 0
#			for line in samfile:
#				if line[0] == '@':
#					continue
#				read_count += 1
#				result = score_read(line, sim_info, analysis_info, args)
#				update_scores(read_scores, result)
#				for sim_match in result.keys():
#					for analysis_match in result[sim_match]:
#						output_writer.writerow(analysis_match)
#	
#	print(f"\nscored {read_count} reads, which were scored: ")
#	pp = pprint.PrettyPrinter(indent = 2)
#	pp.pprint(read_scores)
#
#	write_output_summary(outfile, read_scores, args)

	print(f"saved results to {args.output}")

def process_lines(line_queue, result_queue, sim_info, analysis_info, args):
	"""
	get samfile lines from queue, and write results to file
	"""
	
	while True:
		line = line_queue.get()
		
		# if we're done, stop
		if line is None:
			result_queue.put(None)
			break
		# skip header lines
		elif line[0] == '@':
			continue
		# process this line
		else:
			result = score_read(line, sim_info, analysis_info, args)
			result_queue.put(result)

def read_csv(filename):
	"""
	read the contents of a csv into memory, and return a list of dicts
	where each dict corresponds to a line from the file
	"""
	with open(filename, newline = '') as filehandle:
		
		# create DictReader objects for inputs and read into memory
		reader = csv.DictReader(filehandle, delimiter = '\t')
		data = []
		for row in reader:
			data.append(row)
			
	return data
	
def get_read_properties(line):
	"""
	get properties of read from line of sam file
	"""
	parts = line.split('\t')

	if int(parts[1]) & 64 != 0:
		read_num = "1"
	elif int(parts[1]) & 128 != 0:
		read_num = "2"
	else:
		raise ValueError(f"read {read.qname} is neither read1 nor read2, but reads must be paired")
	
	return {
		'qname' : parts[0],	
		'num' : read_num
	}

def write_output_summary(outfile, read_scores, args):
	"""
	write summary of counts of tp, fp, fn, tn for each type of scoring to file
	"""
	header = ['sim_info_file', 'sim_sam_file', 'analysis_info_file', 'results_file', 'junc_type', 'score_type', 
			  'true_positives', 'true_negatives', 'false_positives', 'false_negatives']
			  
	filenames = [args.sim_info, args.sim_sam, args.analysis_info, args.output]
	types = ['tp', 'tn', 'fp', 'fn']
			  
	with open(args.output_summary, "w") as outfile:
		outfile.write("\t".join(header) + "\n")
		
		for score_type in read_scores:
			for junc_type in read_scores[score_type]:
				if junc_type == 'discord':
					scores = [str(read_scores[score_type][junc_type][type]/2) for type in types]
				else:
					scores = [str(read_scores[score_type][junc_type][type]) for type in types]
				line = filenames + [junc_type, score_type] + scores
				outfile.write("\t".join(line) + "\n")
				
def print_scores(read_scores): 
	
	# recall/sensitvity/true positive rate = (TP)/(TP + FN) - fraction of correct positives
	recall = (read_scores['tp'])/(read_scores['tp'] + read_scores['fn'])
	print(f"true positive rate (TPR/recall/sensitivity)\n  fraction of all positives that are true positives: {recall}")
		
	# specificity/selectivity/true negative rate = (TN)/(TN + FP) - fraction of correct positives
	specificity = (read_scores['tn'])/(read_scores['tn'] + read_scores['fp'])
	print(f"true negative rate (TNR/specificity/selectivity)\n  fraction of all negatives that are true negatives: {specificity}")
		
	# balanced accuracy = (TPR + TNR) / 2	
	print(f"balanced accuracy \n  (TPR + TNR) / 2: {(specificity + recall) / 2}")
		
	# accuracy = (TP + TN) / (TP + TN + FP + FN)
	accuracy =  (read_scores['tn'] + read_scores['tp'] )/ sum([score for type, score in read_scores.items()])
	print(f"accuracy \n  (TP + TN) / (TP + TN + FP + FN): {accuracy}")	
	print()						

def write_results(output_file_header, result_queue, read_scores, n_workers, args):
	"""
	write results from one read to file, and aggregate to read_scores
	"""
	
	# open output file
	outfile =  open(args.output, 'w', newline = '')
	output_writer = csv.DictWriter(outfile, delimiter = '\t', 
										fieldnames = output_file_header)
	output_writer.writeheader()

	none_count = 0
	while True:
		result = result_queue.get()
		# check that there are some results for this read
		if result is None:
			none_count += 1
			if none_count == n_workers:
				pp = pprint.PrettyPrinter(indent = 2)
				pp.pprint(read_scores)
				write_output_summary(output_writer, read_scores, args)
				break
		else:
			# update read_scores and 
			update_scores(read_scores, result)
			for sim_match in result.keys():
				for analysis_match in result[sim_match]:
					output_writer.writerow(analysis_match)
	
	outfile.close()				
	return read_scores

def update_scores(read_scores, result):
	"""
	update scores based on the results from one read, and write the results to file
	"""

	# check that there are some results for this read
	assert len(result) > 0
	# write this read to output
	for sim_match in result.keys():

		# get type (discordant or chimeric)
		junc_type = sim_match.split('_')[2]

		for analysis_match in result[sim_match]:
			
			# get each score type
			for score_type in read_scores:
				score = analysis_match[score_type]
				read_scores[score_type][junc_type][score] += 1			
		
def score_read(line, sim_info, analysis_info, args):
	# score read in terms of positive/negate and true/false
	# simplest scoring is based on whether or not a read was found in the analysis results
	# and whether or not it should have been found (based on simulation results)
	
	read = get_read_properties(line)

	analysis_matches = {}

	# find any simulated integrations involving this read
	# there might be multiple integrations crossed by a read or pair
	read_sim_matches = look_for_read_in_sim(read['qname'], read['num'], sim_info)

	# each read can be involved in one or more chimeric junctions, and at most one 
	# discordant junctions.  if read is not involved in at least one chimeric and one
	# discordant junction, add empty entries to read_sim_matches
	if not any([key.split('_')[2] == 'discord' for key in read_sim_matches.keys()]):
		read_sim_matches['__discord'] = None
	if not any([key.split('_')[2] == 'chimeric' for key in read_sim_matches.keys()]):
		read_sim_matches['__chimeric'] = None

	#  if read_sim_matches is not empty, it crosses at least one integration
	# read is either true positive or false negative
	assert len(read_sim_matches) >= 2
	for int_id, sim_row in read_sim_matches.items():

		# for each integration crossed by a read or pair, it could be in the results more than once
		# as a pair (for one integration) and as a soft-clipped read (for another integration)
		analysis_matches[int_id] = look_for_read_in_analysis(read['qname'], read['num'], int_id, sim_row, analysis_info, args)
		
	score_matches(analysis_matches, args)
	
	return analysis_matches

def score_matches(analysis_matches, args):
	"""
	for each match between analysis and simulation for a read, decide if
	the match is a true positive, false positive, true negative or false negative
	
	assign several scores:
	'found_score' is purely based on read IDs - if the read crosses an integration in the 
	simulation, and if it was found in the analysis results
	
	'host_score' takes into account where the integration was in the host: to be a true
	positive, read must have a 'found_score' of 'tp', and additionally must be on the
	correct chromosome and within args.host_threshold bases of where it should be.  A
	read with a 'found_score' of 'tp' that doesn't meet these additional criteria is 'fn'
	
	'virus_score' is the same as 'host_score' but for the location in the virus
	"""
	
	for match_type in analysis_matches.values():

		for match in match_type:
			# if read crosses an integration
			if match['intID'] is not None:
				# if integration was found, and should have been found
				if match['found'] is True:
					match['found_score'] = 'tp'
					
					# check if integration was found in the correct place in the host
					if correct_pos(match, args, ref = 'host'):
						match['host_score'] = 'tp'
					else:
						match['host_score'] = 'fn'							
						
					# check if integration was found in the correct virus
					if correct_pos(match, args, ref = 'virus'):
						match['virus_score'] = 'tp'
					else:
						match['virus_score'] = 'fn'
					
				else:
					match['found_score'] = 'fn'
					match['host_score'] = 'fn'
					match['virus_score'] = 'fn'
			# if read doesn't cross integration
			else:
				if match['found'] is True:
					match['found_score'] = 'fp'
					match['host_score'] = 'fp'
					match['virus_score'] = 'fp'
				else:
					match['found_score'] = 'tn'
					match['host_score'] = 'tn'
					match['virus_score'] = 'tn'
				
	return analysis_matches

def correct_pos(match, args, ref):
	"""
	return true if match is in the correct position (overlap between sim and analysis coordinates), 
	and false otherwise
	
	optionally, can set a threshold (args.wiggle_room) for checking
	the location of the start and stops relative to each other
	
	if the 'type' is chimeric, and a chimeric threshold is set, require both analysis start and stop to be
	within threshold bases of their sim counterparts
	if the 'type' is discordant, and a chimeric threshold is set, require either analysis start and stop to be
	within threshold bases of their sim counterparts
	"""
	
	assert ref in ('host', 'virus')
	if ref == 'host':
		contig_col = 'correct_host_chr'
		overlap_col = 'host_coords_overlap'
		start_col = 'host_start_dist'
		stop_col = 'host_stop_dist'
	else:
		contig_col = 'correct_virus'
		overlap_col = 'virus_coords_overlap'
		start_col = 'virus_start_dist'
		stop_col = 'virus_stop_dist'		
	
	assert isinstance(match, dict)
	assert 'type' in match
	assert contig_col in match
	assert isinstance(match[contig_col], bool)
	assert match['type'] in ['chimeric', 'discord']
	
	# if the match was on the wrong chromosome, it can't be in the correct place
	if match[contig_col] is False:
		return False
		
	# get threshold for match type (chimeric or discordant)
	if match['type'] == 'chimeric':
		threshold = args.wiggle_room
	else:
		threshold = None

	if isinstance(threshold, int):
		assert threshold >= 0
	else:
		assert threshold is None
	
	# if we don't care about the number of bases distance
	# (we never care for discordant pairs)
	if threshold is None:
		return match[overlap_col]
	
	# are the start and stop positions correct?
	start_correct = (match[start_col] <= threshold)
	stop_correct = (match[stop_col] <= threshold)
	
	# return true if conditions are met for correct position, false otherwise
	if match[contig_col] and start_correct and stop_correct and match[overlap_col]:
		return True

	return False	

def look_for_read_in_sim(readID, read_num, sim_info):
	"""
	given a read id, find any integrations which involved this read
	"""
	
	sim_ints = {}
	
	# look through rows of sim info for matches
	for sim_row in sim_info:
		
		# look in chimeric
		if f"{readID}/{read_num}" in sim_row['left_chimeric'].split(";"):
			sim_ints[f"{sim_row['id']}_left_chimeric"] = sim_row
		
		if f"{readID}/{read_num}" in sim_row['right_chimeric'].split(";"):
			sim_ints[f"{sim_row['id']}_right_chimeric"] = sim_row
			
		# look in discordant
		if readID in sim_row['left_discord'].split(";"):
			sim_ints[f"{sim_row['id']}_left_discord"] = sim_row
			
		if readID in sim_row['right_discord'].split(";"):
			sim_ints[f"{sim_row['id']}_right_discord"] = sim_row
			
	return sim_ints

def look_for_read_in_analysis(readID, read_num, int_descr, sim_row, analysis_info, args):
	"""
	Given a read id from a row in the simulation results, try to find a corresponding row
	in the analysis results
	
	Check whether read is on the same side (left or right) in analysis results and simulated info
	as well as if type (discordant or chimeric) is the same
	
	TODO: incorporate info about whether or not it was at the correct location
	
	"""
	# get id side and type from int_descr
	if int_descr is not None:
		id, side, type = int_descr.split('_')
		if id == '':
			id = None
		if side == '':
			side = None
	else:
		id, side, type = (None, None, None)
	
	# side and type will be None if the read wasn't found in the simulation information
	# but we're looking for it anyway in the analysis results
	assert side in ['left', 'right', None]
	assert type in ['chimeric', 'discord', None]
	
	# if the type is chimeric, we're looking only for a read with the read number appended
	if type == 'chimeric':
		readID = f"{readID}/{read_num}"
		
	# does this read cross an integration?
	cross_int = (sim_row is not None)
	
	# dictionary to store matches for this read
	sim_matches = {'readID' : readID,
					'intID' : id,
					'found' : False,
					'n_found' : 0,
					'side'  : side,
					'correct_side' : None,
					'type'  : type,
					'correct_type' : None,
					'correct_host_chr' : None,
					'host_start_dist' : None, 
					'host_stop_dist' : None, 
					'correct_virus' : None,
					'virus_start_dist' : None,
					'virus_stop_dist' : None, 
					'ambig_diff' : None
					}
	matches = []
	
	
	# look through rows of analysis for matches
	for analysis_row in analysis_info:
	
		# check for  readID
		if analysis_row['ReadID'] == readID:
			
			sim_matches['found'] = True
			sim_matches['n_found'] = 1
			
			# check for correct side
			if analysis_row['Orientation'] == 'hv':
				analysis_side = 'left'
			elif analysis_row['Orientation'] == 'vh':
				analysis_side = 'right'
			else:
				raise ValueError(f'unknown Orientation in analysis results for read {readID}')
			sim_matches['correct_side'] = (analysis_side == side)
			
			# check for correct type
			if analysis_row['OverlapType'] == 'discordant':
				analysis_type = 'discord'
			elif analysis_row['OverlapType'] in ['gap', 'overlap', 'none']:
				analysis_type = 'chimeric'
			else:
				raise ValueError(f'unknown OverlapType in analysis results for read {readID}')
			sim_matches['correct_type'] = (analysis_type == type)
			
			#if this read crosses a simulated int, check for matches between sim and analysis properties
			if cross_int:
				# check for correct host chromosome, 
				sim_matches['correct_host_chr'] = (analysis_row['Chr'] == sim_row['chr'])
			
				# check distance between sim and analysis integration sites in host
				if sim_matches['correct_host_chr']:
					if side == 'left':
						sim_start = int(sim_row['hPos'])
						sim_ambig = int(sim_row['juncLengths'].split(',')[0])
						sim_stop = sim_start + sim_ambig
					else:
						sim_left_ambig = int(sim_row['juncLengths'].split(',')[0])
						sim_right_ambig = int(sim_row['juncLengths'].split(',')[1])
						sim_start = int(sim_row['hPos']) + sim_left_ambig + int(sim_row['hDeleted'])
						sim_stop = sim_start + sim_right_ambig						
						
					analysis_start = int(analysis_row['IntStart'])
					analysis_stop = int(analysis_row['IntStop'])
					
					# check distance between starts and stops	
					sim_matches['host_start_dist'] = abs(sim_start - analysis_start)
					sim_matches['host_stop_dist'] = abs(sim_stop - analysis_stop)
					
					# check if analysis and sim coordates overlap (allowing some wiggle room)
					analysis_start -= args.wiggle_room
					analysis_stop += args.wiggle_room
					sim_matches['host_coords_overlap'] = intersect(sim_start, sim_stop, analysis_start, analysis_stop)	
			
				# check for correct virus
				sim_matches['correct_virus'] = (analysis_row['VirusRef'] == sim_row['virus'])
				
				# check viral coordinates are correct
				if sim_matches['correct_virus']:
					sim_start, sim_stop = get_virus_coordinates(sim_row, side)
					analysis_start = int(analysis_row['VirusStart'])
					analysis_stop = int(analysis_row['VirusStop'])
					
					# check distance between starts and stops
					sim_matches['virus_start_dist'] = abs(sim_start - analysis_start)
					sim_matches['virus_stop_dist'] = abs(sim_stop - analysis_stop)
					
					# check if analysis and sim coordates overlap (allowing some wiggle room)
					analysis_start -= args.wiggle_room
					analysis_stop += args.wiggle_room
					sim_matches['virus_coords_overlap'] = intersect(sim_start, sim_stop, analysis_start, analysis_stop)
					
								
			# append a copy of sim_matches, so that in the case of multiple matches
			# we can check for the best one
			matches.append(dict(sim_matches))

	# if we didn't get any matches
	if len(matches) == 0:
		matches.append(sim_matches)
		
	# if we found more than one match, need to update n_found in each
	if len(matches) > 1:
		print(f"WARNING: read {readID} found more than once in analysis matches")
		for match_dict in matches:
			match_dict['n_found'] = len(matches)
	
	return matches
	
def intersect(start1, stop1, start2, stop2):
	"""
	check if two intervals, defined by start1, stop1, start2 and stop2, intersect
	
	use 0-based coordinates, but allow book-ending, so (3, 4) and (4, 5) do intersect
	return True if they do, otherwise return False
	
	"""
	assert isinstance(start1, int)
	assert isinstance(stop2, int)
	assert isinstance(start2, int)
	assert isinstance(stop2, int)
	assert start1 <= stop1
	assert start2 <= stop2
	
	# if interval 1 is completely to the left of interval 2
	if stop1 < start2:
		return False
	
	# if interval 1 is completely to the right of interval2
	if stop2 < start1:
		return False
		
	return True
	
def read_vifi_file(filename):	
	"""
	read output from ViFi (HpvIntegrationInfo.tsv), and convert into a format similar to our pipeline
	so that it can be used in the same way
	
	since reads are listed on individal lines after the location of the site, collect all the 
	lines with read IDs, and match up the viral and host information
	
	note that ViFi doesn't indicate if the read was read 1 or read 2, so can't check this
	"""
	data = []
	with open(filename) as filehandle:
		
		next(filehandle) # skip header row 
		site = None
		buffer = []
		
		for row in filehandle:
			# new integration site
			if row[:3] == "##=":
				# process previous integration rows
				data += create_vifi_rows(buffer, site)
				
				# get information about this site
				buffer = []
				site = next(filehandle).strip()
				
			# read for current integration site
			else:
				buffer.append(row.strip())
	
	# process last integration
	data += create_vifi_rows(buffer, site)
			
	pdb.set_trace()
			
	return data	

def create_vifi_rows(buffer, site):
	"""
	use lines of ViFi results file to produce dict in our output style, for use in scoring
	site contains the line with information about the integration site location
	buffer contains a buffer of the lines with information about each read
	
	note that ViFi doesn't indicate if the read was read 1 or read 2, so can't check this
	calling all the reads 'discordant', because we don't know which they are
	"""
	# for the first integration, we haven't got any information yet
	if site is None:
		return []
	
	# get information about integration site
	site = site.split('\t')
	chr = site[0]
	data = []
	
	# pair up host and virus information for each read
	pairs = {}
	for row in buffer:
		info = row.split('\t')
		readID = info[0][2:]
		pos = info[2]
		forward = (info[3] == "True")
		ref = info[1]
		viral = info[1] != chr
		
		if readID not in pairs:
			pairs[readID] = {}
		
		if viral:
			key = 'viral'
		else:
			key = 'host'

		pairs[readID][key] = {}
		pairs[readID][key]['ref'] = ref
		pairs[readID][key]['pos'] = pos
		pairs[readID][key]['forward'] = forward
			
	# make one dict for each read
	for readID in pairs:
		read_data = dict()
		read_data['Chr'] = chr
		read_data['IntStart'] = pairs[readID]['host']['pos']
		read_data['IntStop'] = pairs[readID]['host']['pos']
		read_data['VirusRef'] = pairs[readID]['viral']['ref']
		read_data['VirusStart'] =  pairs[readID]['viral']['pos']
		read_data['VirusStop'] = pairs[readID]['viral']['pos']
		read_data['Orientation'] = 'hv' if pairs[readID]['host']['forward'] else 'vh'
		read_data['OverlapType'] = 'none'
		read_data['Type'] = 'discordant'
		read_data['ReadID'] = readID
		data.append(read_data)
	
	return data
	
def read_polyidus_csv(filename):
	"""
	read output from polyidus (HpvIntegrationInfo.tsv), and convert into a format similar to our pipeline
	so that it can be used in the same way
	"""
	with open(filename, newline = '') as filehandle:
		
		# create DictReader objects for inputs and read into memory
		reader = csv.DictReader(filehandle, delimiter = '\t')
		data = []
		read_ids = []
		
		for row in reader:
			row_data = {}
			row_data['Chr'] = row['ChromHost']
			row_data['VirusRef'] = row['ChromViral']
			row_data['OverlapType'] =  'none'
			row_data['Type'] = 'chimeric'
			
			## TODO - each row combines muliple reads, each with a position in the host and virus, and a strand
			# but the number of readnames doesn't match the number of other attributes, so it's unclear
			# how to join up the read names with the other attributes.
			
			hPositions = row['PositionHost'].split(', ')
			vPositions = row['PositionViral'].split(', ')
			hOris = row['StrandHost'].split(', ')
			readIDs = row['ReadNames'].split(', ')
			
			# make one row per read, if we haven't already used this read
			for i, read in enumerate(readIDs):
				if read not in read_ids:
					read_ids.append(read)
					
					
					# need to make copy of dict
					row_data = dict(row_data)
					
					# add info about this read to dict
					row_data['IntStart'] = hPositions[i]
					row_data['IntStop'] = hPositions[i]
					row_data['VirusStart'] = vPositions[i]
					row_data['VirusStop'] = vPositions[i]
					row_data['Orientation'] = 'hv' if hOris[i] == "Positive" else 'vh'
					row_data['type'] = 'chimeric'
					row_data['ReadID'] = read[:-2] + '/' + read[-1]
					
					data.append(row_data)

	return data
	
	
def get_virus_coordinates(sim_row, side):
	"""
	get expected (based on simulation) coordinates of integration in virus/vector
	"""
	
	assert side in ('left', 'right')
	
	# get number of bases at the junction
	if side == 'left':
		n = int(sim_row['juncLengths'].split(",")[0])
	else:
		n = int(sim_row['juncLengths'].split(",")[1])
	
	# convert breakpoints into a usable list of lists
	breakpoints = sim_row['vBreakpoints'].split(";")
	assert all([len(i.split(",")) == 2 for i in breakpoints])
	breakpoints = [(int(i.split(",")[0][1:]), int(i.split(",")[1][:-1])) for i in breakpoints]

	# also need oris
	oris = sim_row['vOris'].split(",")
	assert all([i in( "+", "-") for i in oris])
	
	# get initial start and stop
	if side == 'left':
		i = 0
		piece = breakpoints[i]
		if oris[i] == '+':
			start = piece[0]
			stop = piece[0] + n
		else:
			start = piece[1] - n
			stop = piece[1]
	else:
		i = -1
		piece = breakpoints[i]
		if oris[i] == '+':
			start = piece[1] - n
			stop = piece[1]
		else:
			start = piece[0]
			stop = piece[0] + n	
		
	# check if both start and stop are within the piece - if not print a warning
	# but don't go to the next piece (because the coordinates may have a gap and this
	# creates complexities) 
	
	if start < piece[0] or stop > piece[1]:
		print("warning: start and stop viral coordinates for integration {sim_row['id]} are outside of the first or last piece!  virus_score might be wrong")		
	
	return (start, stop)
	

if __name__ == "__main__":
	t0 = time.time()
	main(argv[1:])
	t1 = time.time()
	print(f"total execution time: {t1-t0}s")
