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

from sys import argv
import argparse
import csv
import pysam
import pdb

def main(argv):
	#get arguments
	parser = argparse.ArgumentParser(description='simulate viral integrations')
	parser.add_argument('--sim-info', help='information from simulation with reads annotated', required = True, type=str)
	parser.add_argument('--analysis-info', help='information from analysis of simulated reads', required = True, type=str)
	parser.add_argument('--sim-bam', help='bam file from ART with simulated reads', required=True, type=str)
	parser.add_argument('--output', help='output file for results', required=False, default='results.tsv')
	args = parser.parse_args()
	
	results = []
	
	# import simulated and pipeline information	
	with open(args.sim_info, newline = '') as sim, open(args.analysis_info, newline='') as analysis, open(args.output, "w", newline = '') as outfile:
		
		# open bam file 
		samfile = pysam.AlignmentFile(args.sim_bam)
		
		# create DictReader objects for inputs
		sim_reader = csv.DictReader(sim, delimiter = '\t')
		analysis_reader = csv.DictReader(analysis, delimiter = '\t')
		
		
		# create DictWriter object for output
		output_writer = csv.DictWriter(outfile, delimiter = '\t', 
										fieldnames = ['readID', 'intID', 'found', 'side', 'type', 'n_found',
										'host_start_dist', 'host_stop_dist', 'virus_start_dist',
										'virus_stop_dist', 'ambig_diff']
										)
		output_writer.writeheader()
		
		# keep count of true and false positives and negatives
		read_scores = {'tp':0, 'tn':0, 'fp':0, 'fn':0}
		read_count = 0
		
		# iterate over reads in samfile and check for corresponding lines in simulation and analysis files
		for read in samfile:
		
			read_count += 1
			# score read in terms of positive/negate and true/false
			# simplest scoring is based on whether or not a read was found in the analysis results
			# and whether or not it should have been found (based on simulation results)
			
			# TODO: incorporate info about whether or not it was at the correct location
		
			# find any simulated integrations involving this read
			if read.is_read1 is True:
				read_num = "1"
			elif read.is_read2 is True:
				read_num = "2"
			else:
				raise ValueError(f"read {read.qname} is neither read1 nor read2, but reads must be paired")
			
			# find out if this read crosses any integrations
			read_sim_matches = look_for_read_in_sim(read.qname, read_num, sim, sim_reader)
			
			#  if read_sim_matches is not empty, it crosses at least one integration
			# read is either true positive or false negative
			if len(read_sim_matches) > 0:
				for int_id, sim_row in read_sim_matches.items():
			
					# for each simulated integration, find read in analysis
					side = int_id.split("_")[1]
					type = int_id.split("_")[2]
					read_analysis_matches = look_for_read_in_analysis(read.qname, sim_row['id'], analysis, analysis_reader, side, type)
			
				# score read as true positive or false negative
				if read_analysis_matches[0]['n_found'] > 0:
					read_scores['tp'] += 1
				else:
					read_scores['fn'] += 1
			
			# if read_sim_matches is empty it doesn't cross any integrations
			# read is either false positve or true negative
			else:
				read_analysis_matches = look_for_read_in_analysis(read.qname, '', analysis, analysis_reader,  None, None)
				# score read ase false positive or true negative
				if read_analysis_matches[0]['n_found'] > 0:
					read_scores['fp'] += 1
				else:
					read_scores['tn'] += 1
					
			for match in read_analysis_matches:
				output_writer.writerow(match)
		
		print(f"scored {read_count} reads, which were scored: {read_scores}")

		# recall/sensitvity/true positive rate = (TP)/(TP + FN) - fraction of correct positives
		recall = (read_scores['tp'])/(read_scores['tp'] + read_scores['fn'])
		print(
		f"true positive rate (TPR/recall/sensitivity)\n\tfraction of all positives that are true positives:\n\t{recall}")
		
		# specificity/selectivity/true negative rate = (TN)/(TN + FP) - fraction of correct positives
		specificity = (read_scores['tn'])/(read_scores['tn'] + read_scores['fp'])
		print(f"true negative rate (TNR/specificity/selectivity)\n\tfraction of all negatives that are true negatives:\n\t{specificity}")
		
		# balanced accuracy = (TPR + TNR) / 2	
		print(f"balanced accuracy \n\t(TPR + TNR) / 2:\n\t{(specificity + recall) / 2}")
		
		# accuracy = (TP + TN) / (TP + TN + FP + FN)
		accuracy =  (read_scores['tn'] + read_scores['tp'] )/ sum([score for type, score in read_scores.items()])
		print(f"balanced accuracy \n\t(TP + TN) / (TP + TN + FP + FN):\n\t{accuracy}")	
		
# 		sim.seek(0)
# 		for sim_row in sim_reader:
# 			# to store the results that we find for each simulated integration
# 			
# 			# find reads that we're looking for in analysis results
# 			left_chimeric = sim_row['left_chimeric'].split(";")
# 			for read in left_chimeric:
# 				results = look_for_read_in_analysis(read, sim_row, analysis, analysis_reader, side = 'left', type='chimeric')
# 				for result in results:
# 					output_writer.writerow(result)
# 				
# 			right_chimeric = sim_row['right_chimeric'].split(";")
# 			for read in right_chimeric:
# 				results = look_for_read_in_analysis(read, sim_row, analysis, analysis_reader, side = 'right', type='chimeric')
# 				for result in results:
# 					output_writer.writerow(result)
# 			
# 			left_discord = sim_row['left_discord'].split(";")
# 			for read in left_discord:
# 				results = look_for_read_in_analysis(read, sim_row, analysis, analysis_reader, side = 'right', type='discord')
# 				for result in results:
# 					output_writer.writerow(result)
# 				
# 			right_discord = sim_row['right_discord'].split(";")
# 			for read in right_discord:
# 				results = look_for_read_in_analysis(read, sim_row, analysis, analysis_reader, side = 'right', type='discord')
# 				for result in results:
# 					output_writer.writerow(result)		
				
	print(f"saved results to {args.output}")

def look_for_read_in_sim(readID, read_num, sim_filehandle, sim_reader):
	"""
	given a read id, find any integrations which involved this read
	"""
	
	sim_ints = {}
	
	# look through rows of sim info for matches
	sim_filehandle.seek(0)
	for sim_row in sim_reader:
		
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

def look_for_read_in_analysis(readID, int_id, analysis_filehandle, analysis_reader, side, type):
	"""
	Given a read id from a row in the simulation results, try to find a corresponding row
	in the analysis results
	
	Check whether read is on the same side (left or right) in analysis reuslts and simulated info
	as well as if type (discordant or chimeric) is the same
	"""
	
	# side and type will be None if the read wasn't found in the simulation information
	# but we're looking for it anyway in the analysis results
	assert side in ['left', 'right', None]
	assert type in ['chimeric', 'discord', None]
	
	# dictionary to store matches for this read
	sim_matches = {'readID' : readID,
					'intID' : int_id,
					'found' : False,
					'n_found' : 0,
					'side'  : '',
					'type'  : ''
					}
	matches = []
	
	# look through rows of analysis for matches
	analysis_filehandle.seek(0)
	for analysis_row in analysis_reader:
	
		# check for the correct readID
		if readID == analysis_row['ReadID']:
				
			sim_matches['found'] = True
			sim_matches['n_found'] = 1
			
			# check for correct side
			if analysis_row['Orientation'] == 'hv':
				analysis_side = 'left'
			elif analysis_row['Orientation'] == 'vh':
				analysis_side = 'left'
			else:
				raise ValueError(f'unknown Orientation in analysis results for read {readID}')
			sim_matches['side'] = (analysis_side == side)
			
			# check for correct type
			if analysis_row['OverlapType'] == 'discordant':
				analysis_type = 'discord'
			elif analysis_row['OverlapType'] in ['gap', 'overlap', 'none']:
				analysis_type = 'chimeric'
			else:
				raise ValueError(f'unknown OverlapType in analysis results for read {readID}')
			sim_matches['type'] = (analysis_type == type)
			
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
	main(argv[1:])
