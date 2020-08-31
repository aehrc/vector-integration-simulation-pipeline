#!/usr/bin/env python3


##### NEED TO CHANGE TO SCORE DISCORDANT AND CHIMERIC SEPARATELY


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


# buffer of lines to print to outfile
line_buffer = []
result_buffer = []
# number of lines to keep in buffer before printing
n_buffer = 1000

def main(argv):
	#get arguments
	parser = argparse.ArgumentParser(description='simulate viral integrations')
	parser.add_argument('--sim-info', help='information from simulation with reads annotated', required = True, type=str)
	parser.add_argument('--analysis-info', help='information from analysis of simulated reads', required = True, type=str)
	parser.add_argument('--sim-sam', help='sam file from ART with simulated reads (must not be binary!)', required=True, type=str)
	parser.add_argument('--chimeric-threshold', help='maximum distance a chimeric read integration can be from where it should be before it is scored as incorrect', type=int, default=5)	
	parser.add_argument('--discordant-threshold', help='maximum distance a discordant read pair integration can be from where it should be before it is scored as incorrect', type=int, default=150)
	parser.add_argument('--output', help='output file for results', required=False, default='results.tsv')
	parser.add_argument('--output-summary', help='file for one-line output summary', required=False, default='results-summary.tsv')
	parser.add_argument('--threads', help='number of threads to use', required=False, type=int)
	
	args = parser.parse_args()
	
	results = []
	
	# read simulated and pipeline information into memory
	print(f"opening simulated information: {args.sim_info}")
	with open(args.sim_info, newline = '') as sim:
		
		# create DictReader objects for inputs and read into memory
		sim_reader = csv.DictReader(sim, delimiter = '\t')
		sim_info = []
		for row in sim_reader:
			sim_info.append(row)
	
	print(f"opening analysis information: {args.analysis_info}")
	with open(args.analysis_info, newline='') as analysis:
		analysis_reader = csv.DictReader(analysis, delimiter = '\t')
		analysis_info = []
		for row in analysis_reader:
			analysis_info.append(row)
		
	print(f"opening sam file: {args.sim_sam}")
	with open(args.output, "w", newline = '') as outfile, open(args.sim_sam) as samfile:
		
		
		# create DictWriter object for output
		score_types = ['found_score', 'host_score', 'virus_score']
		header =  ['readID', 'intID'] + score_types + ['found', 'n_found',
					'side', 'correct_side', 'type', 'correct_type', 'correct_host_chr',
					'host_start_dist', 'host_stop_dist', 'correct_virus', 'virus_start_dist',
					'virus_stop_dist', 'ambig_diff']
		output_writer = csv.DictWriter(outfile, delimiter = '\t', 
										fieldnames = header)
		
		output_writer.writeheader()
		
		# keep count of true and false positives and negatives
		read_scores = {'tp':0, 'tn':0, 'fp':0, 'fn':0}
		read_scores  = {score_type : {'chimeric' : dict(read_scores), 'discord' : dict(read_scores)} for score_type in score_types}
		
		read_count = 0
		if read_count % n_print == 0:
			print(f"processed {read_count} reads")
		
		# create pool of worker processes if we're using more than one thread
		if args.threads is None:
			args.threads = cpu_count()
		if args.threads > 1:
			pool = mp.Pool(processes = args.threads)
		
		# callback function for pool can only take one argument, so use functools.Partial
		# to provide read_scores and outfile handle
		callback_func = functools.partial(get_results, output_writer, read_scores)
		
		# iterate over reads in samfile and check for corresponding lines in simulation and analysis files
		for line in samfile:
		
			# skip header lines
			if line[0] == '@':
				continue
			
			# append line to buffer
			read_count += 1
			line_buffer.append(line)
			
			# if the buffer is full
			if len(line_buffer) > n_buffer:
				
				# score reads in buffer
				if args.threads > 1:
					res = pool.apply_async(score_read, 
						args = (line_buffer, sim_info, analysis_info, args), 
						callback = callback_func
					   )
				else:
					results = score_read(line_buffer, sim_info, analysis_info, args)
					callback_func(results)
				
				# clear buffer
				line_buffer.clear()
				
		
		# score last buffer
		if args.threads > 1:
			res = pool.apply_async(score_read, 
				args = (line_buffer, sim_info, analysis_info, args), 
				callback = callback_func
				)
			pool.close()
			pool.join()
		else:
			results = score_read(line_buffer, sim_info, analysis_info, args)
			callback_func(results)
		
		
	#print(f"\nscored {read_count} reads, which were scored: ")
	pp = pprint.PrettyPrinter(indent = 2)
	pp.pprint(read_scores)
	
	write_output_summary(outfile, read_scores, args)
	print(f"saved results to {args.output}")

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

def get_results(output_writer, read_scores, result_buffer):
	"""
	update scores based on the results from one read, and write the results to file
	"""
	# check that there are some results
	assert len(result_buffer) > 0			

	for result in result_buffer:
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
				# write row
				output_writer.writerow(analysis_match)
	result_buffer.clear()
		
def score_read(line_buffer, sim_info, analysis_info, args):
	# score read in terms of positive/negate and true/false
	# simplest scoring is based on whether or not a read was found in the analysis results
	# and whether or not it should have been found (based on simulation results)

	for line in line_buffer:
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
			analysis_matches[int_id] = look_for_read_in_analysis(read['qname'], read['num'], int_id, sim_row, analysis_info)
				
		# score matches
		result_buffer.append(score_matches(analysis_matches, args))
	
	return result_buffer	

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
					start_correct = (match['host_start_dist'] <= args.chimeric_threshold)
					stop_correct = (match['host_stop_dist'] <= args.chimeric_threshold)
					if match['correct_host_chr'] and start_correct and stop_correct:
						match['host_score'] = 'tp'
					else:
						match['host_score'] = 'fn'
						
					# check if integration was found in the correct place in the host
					if match['correct_virus'] is True:
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

def look_for_read_in_analysis(readID, read_num, int_descr, sim_row, analysis_info):
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
				analysis_side = 'left'
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
						sim_start = int(sim_row['hPos']) + sim_left_ambig
						sim_stop = sim_start + sim_right_ambig						
						
					analysis_start = int(analysis_row['IntStart'])
					analysis_stop = int(analysis_row['IntStop'])
						
					sim_matches['host_start_dist'] = abs(sim_start - analysis_start)
					sim_matches['host_stop_dist'] = abs(sim_stop - analysis_stop)	
			
				# check for correct virus
				sim_matches['correct_virus'] = (analysis_row['VirusRef'] == sim_row['virus'])
			
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

if __name__ == "__main__":
	t0 = time.time()
	main(argv[1:])
	t1 = time.time()
	print(f"total execution time: {t1-t0}s")
