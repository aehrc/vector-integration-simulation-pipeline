#!/usr/bin/env python3

# Score simulated vs detected integrations on a per-integration basis

# do this by counting the number of lines in found-info file that 
# fall within a specified window around each integration in sim-info file

# integrations may have deletions from host and ambiguous bases at each junction,
# so a full consideration of these factors would require scoring each junction (host/virus and virus/host)
# for each integration. however, not all tools indicate which junction, and some (eg vifi) combine both junctions
# in the same 'integration site'.  So instead attempt to consider orientiation, by using different coordinates 
# for the left and right junctions.  If the tool does not output orientation information, fall back to just the 
# position of the integration (without considering ambiguous bases and deletions from the host)

from sys import argv
import argparse
import csv
import pdb

# currently supported tools
supported_tools = ('pipeline', 'vifi', 'polyidus', 'verse', 'seeksv')
score_types = ('overlap', 'coords_min', 'coords_mean')
output_fieldnames = ('id', 'score', 'type', 'chr', 'pos', 'start', 'stop', 'host_del', 
										'dist', 'read_count', 'reads',
										'left_start', 'left_stop', 'left_dist', 'left_reads', 'left_read_count', 
										'right_start', 'right_stop', 'right_dist','right_reads', 'right_read_count')

def main(argv):
	#get arguments
	parser = argparse.ArgumentParser(description='simulate viral integrations')
	parser.add_argument('--sim-info', help='information from simulation', required = True, type=str)
	parser.add_argument('--found-info', help='information from found of simulated reads', required = True, type=str)
	parser.add_argument('--window', help='look this amount either side of simulated hPos for integrations in found', required=False, type=int, default=100)
	parser.add_argument('--output', help='output file', required=False, default='results.tsv')
	parser.add_argument('--output-sim', help='output all simulated integrations with distance to nearest result', required=False)
	parser.add_argument('--output-found', help='output all results with distance to nearest simulated integration', required=False)	
	parser.add_argument('--summary', help='output summary file', required=False, default='summary.tsv')	
	parser.add_argument('--analysis-tool', choices=supported_tools, required=True)
	parser.add_argument('--score-type', choices=score_types, required=False, default='coords_mean')
	
	args = parser.parse_args()
	
	# read in information about which integrations were found
	if args.analysis_tool == "pipeline":
		found = parse_results_tsv(args.found_info)
	elif args.analysis_tool == "polyidus":
		found = parse_polyidus(args.found_info)
	elif args.analysis_tool == "vifi":
		found = parse_vifi(args.found_info)
	elif args.analysis_tool == "verse":
		found = parse_verse(args.found_info)
	elif args.analysis_tool == "seeksv":
		found = parse_seeksv(args.found_info)
	
	# read in infomration about simulated integrations
	sim = parse_tsv(args.sim_info)
	
	# remove any simulated integrations that don't have any supporting reads
	sim = remove_unsupported_ints(sim)
	
	# to store results
	scored_sim_results = []
	
	# to find true positives and false negatives, loop over rows of sim
	for sim_int in sim:
		
		# get number of ambiguous bases at each junction
		left_ambig = int(sim_int['juncLengths'].split(",")[0])
		right_ambig = int(sim_int['juncLengths'].split(",")[1])
		host_del = int(sim_int['hDeleted'])
		sim_pos = int(sim_int['hPos'])
		
		# to store the results that we find for each simulated integration

		sim_result = {'id' : sim_int['id'],
									'chr': sim_int['chr'],
									'pos': sim_pos,
									'host_del': host_del,
									'start': sim_pos,
									'stop' : sim_pos,
									'left_start': sim_pos,
									'left_stop' : sim_pos + left_ambig,
									'right_start': sim_pos + host_del + left_ambig, 
									'right_stop' : sim_pos + host_del + left_ambig + right_ambig, 
									'dist' : None, 
									'left_dist': None,
									'right_dist': None,
									'left_read_count': 0,
									'right_read_count': 0,
									'read_count': 0,
									'left_reads': [],
									'right_reads': [],
									'reads': [],
									'score' : '',
									'type' : 'sim'
									}

		sim_result = find_sim_in_found(sim_result, args, found)
		scored_sim_results.append(sim_result)

	# write these results to file, if desired
	if args.output_sim is not None:
		with open(args.output_sim, "w", newline = '') as outfile:
			# create DictWriter object for scored_sim_results
			output_writer = csv.DictWriter(outfile, delimiter = '\t', 
											fieldnames = output_fieldnames)
			output_writer.writeheader()
			for row in scored_sim_results:
				output_writer.writerow(row)

	# check found for integrations that we saw but didn't expect to
	found_results = []
	for result in found:
		# find closest simulated integration to this output integration
		sim_found_result = check_result_in_sim(result, scored_sim_results, args)
		
		# if this result wasn't in the sim integrations
		if 'found' not in result:
	
			scores['fp'] += 1
			sim_found_result['score'] = 'fp'
			found_results.append(sim_found_result)
			
		else:
			sim_found_result['score'] = 'tp'
			found_results.append(sim_found_result)

	# output these results, if user asked for them
	if args.output_found is not None:
		with open(args.output_found, "w", newline = '') as outfile:
			# create DictWriter object for scored_sim_results
			output_writer = csv.DictWriter(outfile, delimiter = '\t', 
											fieldnames = output_fieldnames)
			output_writer.writeheader()
			for row in found_results:
				output_writer.writerow(row)	

	# only add the false positives to the combined output
	scored_sim_results += [res for res in found_results if res['score'] == 'fp']
			
	# write results to file
	with open(args.output, "w", newline = '') as outfile:
		# create DictWriter object for scored_sim_results
		output_writer = csv.DictWriter(outfile, delimiter = '\t', 
										fieldnames = output_fieldnames)
		output_writer.writeheader()
		for row in scored_sim_results:
			row['reads'] = ";".join(row['reads'])
			row['left_reads'] = ";".join(row['left_reads'])
			row['right_reads'] = ";".join(row['right_reads'])
			output_writer.writerow(row)

	
	# check scored_sim_results (which has one entry for each simulated integration) 
	# for true positives and false negatives	
	scores = {'sim_info'      : args.sim_info,
			  		'found_info' : args.found_info,
			  		'analysis_tool' : args.analysis_tool,
			  		'window'		  : args.window,
			  		'coords_score_type': args.score_type,
			  		'tp':0, 'fp':0, 'tn':0, 'fn':0 }

	# add scores up
	for result in scored_sim_results:
		# integrations that we saw, and expected to see
		for score_type in ('tp', 'fn', 'fp'):
			if result['score'] == score_type:
				scores[score_type] += 1

	
	# write summary
	with open(args.summary, "w", newline = '') as outfile:
		# create dictwrite object for summary
		output_writer = csv.DictWriter(outfile, delimiter = '\t',
										fieldnames = ('sim_info', 'found_info',
													   'analysis_tool', 'window', 'coords_score_type',
													   'tp', 'tn', 'fp', 'fn'))
		output_writer.writeheader()
		output_writer.writerow(scores)
		
	print(f"true positives: {scores['tp']}")	
	print(f"true negatives: {scores['tn']}")
	print(f"false positives: {scores['fp']}")	
	print(f"false negatives: {scores['fn']}")
	print(f"\nsaved output to {args.output}\nsaved summary to {args.summary}")
	

def check_result_in_sim(result, sim, args):

		closest_dist, closest = get_closest_dist(result, sim, args)

		result_dict = {'id' : result['id'],
									'chr': result['Chr'],
									'pos': f"{result['IntStart']}/{result['IntStop']}",
									'host_del': None,
									'start': result['IntStart'],
									'stop' : result['IntStop'],
									'left_start': None,
									'left_stop' : None,
									'right_start': None, 
									'right_stop' : None, 
									'dist' : closest_dist, 
									'left_dist': None,
									'right_dist': None,
									'left_read_count': 0,
									'right_read_count': 0,
									'read_count': len(result['ReadID']) if len(result['ReadID']) > 0 else 1,
									'left_reads': None,
									'right_reads': None,
									'reads': ";".join(result['ReadID']),
									'type' : 'result'
									}
		return result_dict

def get_closest_dist(result, sim, args):
	"""
	get the simulated integration from list sim that is closest to the result
	and also get the distance between the simulated integration and the result
	"""
	d = 1e10
	for sim_int in sim:
			## check for different chromosomes
		if result['Chr'] != sim_int['chr']:
			continue
		
		# get distance between intervals
		left_d = distance_between_intervals(int(result['IntStart']), int(result['IntStop']), 
																					sim_int['left_start'], sim_int['left_stop'], args.score_type)
		right_d = distance_between_intervals(int(result['IntStart']), int(result['IntStop']), 
																					sim_int['right_start'], sim_int['right_stop'], args.score_type)
		# check if this is closer than before
		if left_d < d:
			d = left_d
			closest = sim_int
		if right_d < d:
			d = right_d
			closest = sim_int
	
	# check if everything was on the wrong chromosome
	if d == 1e10:
		return None, None
	
	return d, sim_int


def find_sim_in_found(sim_result, args, found):
		# look for integrations in found that overlap this window
		
		# to keep track of the closest result to this integration
		closest = {'found': {}, 'side': "", 'dist': 1e10}
		
		for result in found:
			
			# check for correct chromosome
			if result['Chr'] != sim_result['chr']:
				continue
			
			# check if this integration is within specified distance  
			found_start = int(result['IntStart'])
			found_stop = int(result['IntStop'])
			
			left_dist = distance_between_intervals(sim_result['left_start'], sim_result['left_stop'], 
																							found_start, found_stop, args.score_type)
			if left_dist < closest['dist']:
				closest['side'] = 'left'
				closest['found'] = result
				closest['dist'] = left_dist
			left_overlap = (left_dist <= args.window)
			right_dist = distance_between_intervals(sim_result['right_start'], sim_result['right_stop'], 
																							found_start, found_stop, args.score_type)
			if right_dist < closest['dist']:
				closest['side'] = 'right'
				closest['found'] = result	
				closest['dist'] = right_dist
			left_overlap = (left_dist <= args.window)
			right_overlap = (right_dist <= args.window)
			
		
			if left_overlap or right_overlap:
				
				sim_result['score'] = 'tp'
				
				# increment counters
				sim_result['read_count'] += len(result['ReadID']) if len(result['ReadID']) > 0 else 1
				sim_result['reads'] += result['ReadID']
				
				# mark this integration in found
				if 'found' not in result:
					result['found'] = 1
				else:
					result['found'] += 1
					
				# find distance between expected and found start and stop
				if left_overlap:
					if sim_result['left_dist'] is None:
						sim_result['left_dist'] = left_dist
					else:
						sim_result['left_dist'] = min(left_dist, sim_result['left_dist'])
					
					sim_result['left_reads'] += result['ReadID']
				
				if right_overlap:
					if sim_result['right_dist'] is None:
						sim_result['right_dist'] = right_dist
					else:
						sim_result['right_dist'] = min(right_dist, sim_result['right_dist'])
					sim_result['right_reads'] += result['ReadID']
		
		# if this simulated integration wasn't found
		# then we call it a 'false negative'
		if sim_result['read_count'] == 0:
			sim_result['score'] = 'fn'
			sim_result['dist'] = closest['dist'] if closest['dist'] != 1e10 else None
			sim_result['reads'] += closest['found']['ReadID'] if closest['dist'] != 1e10 else None
			
		# otherwise, assign a minimum distance 
		else:
			if sim_result['left_dist'] is None:
				sim_result['dist'] = sim_result['right_dist']
			elif sim_result['right_dist'] is None:
				sim_result['dist'] = sim_result['left_dist']	
			else:
				sim_result['dist'] = min(sim_result['left_dist'], sim_result['right_dist'])
			
		
		sim_result['reads'] = ";".join(sim_result['reads'])
		sim_result['left_reads'] = ";".join(sim_result['left_reads'])
		sim_result['right_reads'] = ";".join(sim_result['right_reads'])
		
		return sim_result	


def distance_between_intervals(start1, stop1, start2, stop2, dist_type):
	"""
	get the distance between two intervals
	if dist_type == 'overlap', the distance is the smallest distance between the two intervals, or 0 if they overlap (bedtools style)
	if dist_type == 'coords', the distance is the smaller of the difference between the two starts and the two stops
	"""
	assert stop1 >= start1
	assert stop2 >= start2
	assert dist_type in score_types
	if dist_type == 'overlap':
		# if interval 1 to the left of interval 2
		if stop1 < start2:
			return start2 - stop1
		# if interaval 1 to the right of interval 2
		if start1 > stop2:
			return start1 - stop2
		# otherwise, overlapping
		return 0
	elif dist_type == 'coords_min':
		d1 = abs(start1 - start2)
		d2 = abs(stop1 - stop2)
		return(min(d1, d2))
	elif dist_type == 'coords_mean':
		d1 = abs(start1 - start2)
		d2 = abs(stop1 - stop2)
		return((d1 + d2) / 2)	
	

def get_start_stop_dist(sim_start, sim_stop, found_start, found_stop, start_dist, stop_dist):
	"""
	get the distance between sim_start and found_start, and between sim_stop and found_stop
	"""
	if start_dist is None:
		start_dist = abs(sim_start - found_start)
		stop_dist = abs(sim_stop - found_stop)
	else:
		start_dist = min(abs(sim_start - found_start), start_dist)
		stop_dist = min(abs(sim_stop - found_stop), stop_dist)
	return start_dist, stop_dist
	
def parse_tsv(path):
	"""
	parse a tsv, and store each row as a dict in a list
	"""
	with open(path, newline = '') as handle:
		reader = csv.DictReader(handle, delimiter = '\t')
		data = []
		
		# get info from rows
		for row in reader:
			
			data.append(row)
			
	return data	

def overlap(a1, a2, b1, b2):
	"""
	return true if the interval that starts at b1 and stops at b2 ovleraps with the interval
	that starts at a1 and ends at a2
	"""
	assert a1 <= a2
	assert b1 <= b2
	assert isinstance(a1, int) and isinstance(a2, int) and isinstance(b1, int) and isinstance(b2, int)
	
	# if a interval is completely to the left of the b interval
	if a2 < b1:
		return False
	# if a interval is completely to the right of the b interval
	elif a1 > b2:
		return False
	else:
		return True
		
def remove_unsupported_ints(sim):
	"""
	remove any integrations that are not supported by any reads
	"""

	for i, row in enumerate(sim):
		if row['left_chimeric'] != "":
			continue
		if row['right_chimeric'] != '':
			continue
		if row['left_discord'] != '':
			continue
		if row['right_discord'] != '':
			continue
		sim.pop(i)
	
	return sim
	
	
def parse_polyidus(path):
	"""
	parse polyidus output file (exactHpvIntegrations.tsv), and produce a data structure 
	similar to parse_merged_bed
	
	assign unique id to each integration, which consists of row number in output file
	"""
	with open(path, newline = '') as found:
		found_reader = csv.DictReader(found, delimiter = '\t')
		data = []
		
		# get info from rows
		row_num = 0
		for row in found_reader:
			
			fragments = row['FragmentName'].split(", ")
			
			result = {
				'Chr' 		  : row['Chrom'],
				'IntStart' 	  : int(row['IntegrationSite']),
				'IntStop' 	  : int(row['IntegrationSite']),
				'Orientation' : 'unknown',
				'ReadID'      : fragments,
				'id'		  : row_num
				}
			data.append(result)
			row_num += 1
			
	return data	
	
def parse_results_tsv(path):
	"""
	parse found-info from our pipeline.  Get information necessary to assess true/false postive/negative integration
	
	assign unique ID to each integration, which consists, of Chr/IntStart/Virus/VirusStart/ReadID
	"""
	with open(path, newline = '') as found:
		found_reader = csv.DictReader(found, delimiter = '\t')
		data = []
		
		# get info from rows
		for row in found_reader:

			# only need to keep fields 'Chr', 'IntStart', 'IntStop', 'Orientation', and 'ReadID'
			info = {'Chr' 			: row['Chr'],
					'IntStart' 		: int(row['IntStart']),
					'IntStop' 		: int(row['IntStop']),
					'Orientation'	: 'unknown',
					'ReadID'     	: row['ReadIDs'].split(','),
					'id'			: row['SiteID']
					}

			data.append(info)
			
	return data
	
def parse_vifi(path):
	"""
	parse vifi output file (exactHpvIntegrations.tsv), and produce a data structure 
	similar to parse_results_tsv
	
		assign unique id to each integration, which consists of the first line in the output for this site (before the reads)
	"""	
	data = []

	with open(path) as handle:
		handle.readline()
		reads = set()
		n_ints = 0
		for line in handle:
			# ignore divider lines
			if line[:5] == "##===":
				continue
			
			# if this is a new integration, collect info and then clear buffer
			elif line[0] != "#":

				# if there are some reads for the last integration
				if len(reads) :
					result = {
										'Chr' 		  	: chrom,
										'IntStart' 	  : start,
										'IntStop' 	  : stop,
										'Orientation' : 'unknown',
										'ReadID'			: reads,
										'id'					: n_ints
						}
					
					data.append(result)		
					n_ints += 1							
				
				# get info about this integration
				
				info = line.split()
				chrom = info[0]
				start = int(info[1])
				stop = int(info[2])
				reads = set()
				
			# otherwise, collect reads
			else:
				# get info from this line
				info = line.split()
				# right now we only care about host reads
				if info[1] != chrom:
					continue				
					
				# get if read is at hv (left) or vh (right) junction
				reads.add(info[0][2:])
				
	# add last integration
	if len(reads) :
		result = {
							'Chr' 		  	: chrom,
							'IntStart' 	  : start,
							'IntStop' 	  : stop,
							'Orientation' : 'unknown',
							'ReadID'			: reads,
							'id'					: n_ints
						}
					
		data.append(result)		
					
	return data	

def parse_verse(path):
	"""
	parse verse output file (integration-sites.txt), and produce a data structure 
	similar to parse_merged_bed
	
	assign unique id to each integration, which consists of row number in output file
	"""
	with open(path, newline = '') as found:
		
		# if file is empty
		if found.readline() == "":
			return []
		
		found.seek(0)
		found_reader = csv.DictReader(found, delimiter = '\t')
		data = []
		
		# get info from rows
		row_num = 0
		for row in found_reader:
			
			result = {
				'Chr' 		  : row['Chromosome 1'] if row['Chromosome 2'] == 'chrVirus' else row['Chromosome 2'],
				'IntStart' 	  : int(row['Position 1'] if row['Chromosome 2'] == 'chrVirus' else row['Position 2']),
				'IntStop' 	  : int(row['Position 1'] if row['Chromosome 2'] == 'chrVirus' else row['Position 2']),
				'Orientation' : 'unknown',
				'ReadID'      : '',
				'id'		  : row_num
				}
			data.append(result)
			row_num += 1
			
	return data	

def parse_seeksv(path):
	"""
	parse seeksv output file, and produce a data structure 
	similar to parse_results_tsv
	
	assign unique id to each integration, which consists of row number in output file
	"""
	with open(path, newline = '') as found:
		
		# if file is empty
		if found.readline() == "":
			return []
		
		found.seek(0)
		found_reader = csv.DictReader(found, delimiter = '\t')
		data = []
		
		# get info from rows
		row_num = 0
		for row in found_reader:
			
			result = {
				'Chr' 		  : row['@left_chr'] if 'chr' in row['@left_chr'] else row['right_chr'] ,
				'IntStart' 	  : int(row['left_pos'] if 'chr' in row['@left_chr'] else row['right_pos']),
				'IntStop' 	  : int(row['left_pos'] if 'chr' in row['@left_chr'] else row['right_pos']),
				'Orientation' : 'unknown',
				'ReadID'      : '',
				'id'		  : row_num
				}
			data.append(result)
			row_num += 1
			
	return data	


if __name__ == "__main__":
	main(argv[1:])
