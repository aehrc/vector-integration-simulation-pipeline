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

from sys import argv
import argparse
import csv
import pysam
import pdb
import time

def main(argv):
	#get arguments
	parser = argparse.ArgumentParser(description='simulate viral integrations')
	parser.add_argument('--sim-info', help='information from simulation with reads annotated', required = True, type=str)
	parser.add_argument('--analysis-info', help='information from analysis of simulated reads', required = True, type=str)
	parser.add_argument('--sim-bam', help='bam file from ART with simulated reads', required=True, type=str)
	parser.add_argument('--output', help='output file for results', required=False, default='results.tsv')
	parser.add_argument('--output-summary', help='file for one-line output summary', required=False, default='results-summary.tsv')

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
		
	with open(args.output, "w", newline = '') as outfile:
		
		# open bam file 
		print(f"opening sam file: {args.sim_bam}")
		samfile = pysam.AlignmentFile(args.sim_bam)
		
		# create DictWriter object for output
		output_writer = csv.DictWriter(outfile, delimiter = '\t', 
										fieldnames = ['readID', 'intID', 'score', 'found', 'n_found',
										'side', 'correct_side', 'type', 'correct_type', 'correct_host_chr',
										'host_start_dist', 'host_stop_dist', 'correct_virus', 'virus_start_dist',
										'virus_stop_dist', 'ambig_diff']
										)

		
		output_writer.writeheader()
		
		# keep count of true and false positives and negatives
		read_scores = {'tp':0, 'tn':0, 'fp':0, 'fn':0}
		read_count = 0
		
		# iterate over reads in samfile and check for corresponding lines in simulation and analysis files
		for read in samfile:
		
			read_count += 1
			
			# score read
			read_analysis_matches = score_read(read, sim_info, analysis_info)
			
			# get results
			read_scores = get_results(read_analysis_matches, output_writer, read_scores)
			
	print(f"scored {read_count} reads, which were scored: {read_scores}")

	# recall/sensitvity/true positive rate = (TP)/(TP + FN) - fraction of correct positives
	recall = (read_scores['tp'])/(read_scores['tp'] + read_scores['fn'])
	print(f"true positive rate (TPR/recall/sensitivity)\n\tfraction of all positives that are true positives:\n\t{recall}")
		
	# specificity/selectivity/true negative rate = (TN)/(TN + FP) - fraction of correct positives
	specificity = (read_scores['tn'])/(read_scores['tn'] + read_scores['fp'])
	print(f"true negative rate (TNR/specificity/selectivity)\n\tfraction of all negatives that are true negatives:\n\t{specificity}")
		
	# balanced accuracy = (TPR + TNR) / 2	
	print(f"balanced accuracy \n\t(TPR + TNR) / 2:\n\t{(specificity + recall) / 2}")
		
	# accuracy = (TP + TN) / (TP + TN + FP + FN)
	accuracy =  (read_scores['tn'] + read_scores['tp'] )/ sum([score for type, score in read_scores.items()])
	print(f"accuracy \n\t(TP + TN) / (TP + TN + FP + FN):\n\t{accuracy}")	
						
	print(f"saved results to {args.output}")
	
	# write summary of tp, tn, fp, fn
	with open(args.output_summary, "w") as outfile:
		outfile.write("\t".join(['sim_info_file', 'sim_bam_file', 'analysis_info_file', 'results_file', 
														'true_positives', 'true_negatives', 'false_positives', 'false_negatives']) + '\n')
		outfile.write('\t'.join([args.sim_info, args.sim_bam, args.analysis_info, args.output, 
															str(read_scores['tp']), str(read_scores['tn']), str(read_scores['fp']), str(read_scores['fn'])]) + '\n')

def get_results(read_analysis_matches, output_writer, read_scores):
	"""
	update scores based on the results from one read, and write the results to an output file
	"""
	# we should always get a result
	assert len(read_analysis_matches) > 0			
	
	
	# write this read to output
	for sim_match in read_analysis_matches.keys():
	
		for analysis_match in read_analysis_matches[sim_match]:
			score = analysis_match['score']
			read_scores[score] += 1
			output_writer.writerow(analysis_match)
	
	return read_scores
	

def score_read(read, sim_info, analysis_info):
	# score read in terms of positive/negate and true/false
	# simplest scoring is based on whether or not a read was found in the analysis results
	# and whether or not it should have been found (based on simulation results)
			
	

	analysis_matches = {}
	
	if read.is_read1 is True:
		read_num = "1"
	elif read.is_read2 is True:
		read_num = "2"
	else:
		raise ValueError(f"read {read.qname} is neither read1 nor read2, but reads must be paired")
	
	# find any simulated integrations involving this read
	# there might be multiple integrations crossed by a read or pair
	read_sim_matches = look_for_read_in_sim(read.qname, read_num, sim_info)
	
	#  if read_sim_matches is not empty, it crosses at least one integration
	# read is either true positive or false negative
	if len(read_sim_matches) > 0:
		for int_id, sim_row in read_sim_matches.items():

			# for each integration crossed by a read or pair, it could be in the results more than once
			# as a pair (for one integration) and as a soft-clipped read (for another integration)
			analysis_matches[int_id] = look_for_read_in_analysis(read.qname, read_num, int_id, sim_row, analysis_info)

	# if read_sim_matches is empty it doesn't cross any integrations
	# read is either false positve or true negative
	else:
		analysis_matches['no_int'] = look_for_read_in_analysis(read.qname, read_num, '', {}, analysis_info)
		
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
			sim_ints[f"{sim_row['id']}_left_chimeric"] = dict(sim_row)
		
		if f"{readID}/{read_num}" in sim_row['right_chimeric'].split(";"):
			sim_ints[f"{sim_row['id']}_right_chimeric"] = dict(sim_row)
			
		# look in discordant
		if readID in sim_row['left_discord'].split(";"):
			sim_ints[f"{sim_row['id']}_left_discord"] = dict(sim_row)
			
		if readID in sim_row['right_discord'].split(";"):
			sim_ints[f"{sim_row['id']}_right_discord"] = dict(sim_row)
	
			
	return sim_ints

def look_for_read_in_analysis(readID, read_num, int_descr, sim_row, analysis_info):
	"""
	Given a read id from a row in the simulation results, try to find a corresponding row
	in the analysis results
	
	Check whether read is on the same side (left or right) in analysis results and simulated info
	as well as if type (discordant or chimeric) is the same
	
	TODO - for each integration crossed by a read or pair, it could be in the results more than once
	as a pair (for one integration) and as a soft-clipped read (for another integration)
	so need to check which integration read is associated with in sim and results when scoring
	
	TODO: incorporate info about whether or not it was at the correct location
	
	"""
	# get id side and type from int_descr
	if int_descr != "":
		id, side, type = int_descr.split('_')
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
	cross_int = (len(sim_row) == 0)
	
	# dictionary to store matches for this read
	sim_matches = {'readID' : readID,
					'intID' : id,
					'found' : False,
					'n_found' : 0,
					'side'  : side,
					'correct_side' : None,
					'type'  : type,
					'correct_type' : None,
					'score' : '',
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
	
		# check for the correct readID
		
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
				sim_matches['correct_host_chr'] = (analysis_row['chr'] == sim_row['chr'])
			
				# check for correct virus
				sim_matches['correct_virus'] = (analysis_row['VirusRef'] == sim_row['virus'])
			
			
			# score as true positive, true negative, false positive or false negative
			if cross_int is True:
				sim_matches['score'] = 'tp'
			else:
				sim_matches['score'] = 'fp'
			
			# append a copy of sim_matches, so that in the case of multiple matches
			# we can check for the best one
			matches.append(dict(sim_matches))

	
	# if we didn't get any matches
	if len(matches) == 0:
		if cross_int is True:
			sim_matches['score'] = 'fn'
		else:
			sim_matches['score'] = 'tn'
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
