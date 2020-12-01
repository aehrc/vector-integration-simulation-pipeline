#!/usr/bin/env python3

# Score simulated vs detected integrations on a per-integration basis

# do this by counting the number of lines in found-info file that 
# fall within a specified window around each integration in sim-info file

from sys import argv
import argparse
import csv
import pdb

# currently supported tools
supported_tools = ('pipeline', 'vifi', 'polyidus', 'verse', 'seeksv')
score_types = ('overlap', 'coords')

def main(argv):
	#get arguments
	parser = argparse.ArgumentParser(description='simulate viral integrations')
	parser.add_argument('--sim-info', help='information from simulation', required = True, type=str)
	parser.add_argument('--found-info', help='information from found of simulated reads', required = True, type=str)
	parser.add_argument('--window', help='look this amount either side of simulated hPos for integrations in found', required=False, type=int, default=100)
	parser.add_argument('--output', help='output file', required=False, default='results.tsv')
	parser.add_argument('--summary', help='output summary file', required=False, default='summary.tsv')	
	parser.add_argument('--analysis-tool', choices=supported_tools, required=True)
	parser.add_argument('--score-type', choices=score_types, required=False, default='coords')
	
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
		
		# to store the results that we find for each simulated integration
		sim_results = {'id' : sim_int['id'],
										'chr': sim_int['chr'],
										'pos': sim_int['hPos'],
										'read_count': 0,
										'reads': []
									}
						
		# get coordinates in host in which to look for found integrations
		sim_start = int(sim_int['hPos']) - args.window
		sim_stop = int(sim_int['hPos']) + args.window
		
		# look for integrations in found that overlap this window
		for result in found:
			
			# check for correct chromosome
			if result['Chr'] != sim_int['chr']:
				continue
			
			# check for overlap
			found_start = int(result['IntStart'])
			found_stop = int(result['IntStop'])
			
			# check for overlap

			has_overlap = False
			if args.score_type == 'overlap':
				if overlap(found_start, found_stop, sim_start, sim_stop):
					has_overlap = True
			elif args.score_type == 'coords':
				if sim_start <= found_start and sim_stop >= found_start:
					if sim_start <= found_stop and sim_stop >= found_stop:
						has_overlap = True
						
			if has_overlap:
				# increment counters
				sim_results['read_count'] += len(result['ReadID']) if len(result['ReadID']) > 0 else 1
				sim_results['reads'] += result['ReadID']
				
				# mark this integration in found
				if 'found' not in result:
					result['found'] = 1
				else:
					result['found'] += 1
		
		# add to scored
		scored_sim_results.append(sim_results)

	# check scored_sim_results (which has one entry for each simulated integration) 
	# for true positives and false negatives	
	scores = {'sim_info'      : args.sim_info,
			  		'found_info' : args.found_info,
			  		'analysis_tool' : args.analysis_tool,
			  		'window'		  : args.window,
			  		'coords_score_type': args.score_type,
			  		'tp':0, 'fp':0, 'tn':0, 'fn':0 }

	for result in scored_sim_results:
		# integrations that we saw, and expected to see
		if result['read_count'] > 0:
			scores['tp'] += 1
			result['score'] = 'tp'
		# integrations that we didn't
		else:
			scores['fn'] += 1
			result['score'] = 'fn'
	
	
	# check found for integrations that we saw but didn't expect to
	for result in found:

		# if this result wasn't in the sim integrations
		if 'found' not in result:
		
			# if this is the first read we've seen from this integration
			scores['fp'] += 1
				
			# each wrong integration is one false positive
			result_dict = {
				'id': result['id'],
				'score': 'fp',
				'chr': result['Chr'],
				'pos': f"{result['IntStart']}/{result['IntStop']}",
				'read_count': len(result['ReadID']) if len(result['ReadID']) > 0 else 1,
				'reads': result['ReadID']
			}
			scores['fp'] += 1
			scored_sim_results.append(result_dict)
					
				
	# join read lists for output
	for result_dict in scored_sim_results:
		result_dict['reads'] = ";".join(result_dict['reads'])
			
	# write results to file
	with open(args.output, "w", newline = '') as outfile:
		# create DictWriter object for scored_sim_results
		output_writer = csv.DictWriter(outfile, delimiter = '\t', 
										fieldnames = ('id', 'score', 'chr', 'pos',
										'read_count', 'reads'))
		output_writer.writeheader()
		for row in scored_sim_results:
			output_writer.writerow(row)

	
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
			
			fragments = ";".join(row['FragmentName'].split(", "))
			
			result = {
				'Chr' 		  : row['Chrom'],
				'IntStart' 	  : row['IntegrationSite'],
				'IntStop' 	  : row['IntegrationSite'],
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
					'IntStart' 		: row['IntStart'],
					'IntStop' 		: row['IntStop'],
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
										'ReadID'			: ";".join(reads),
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
							'ReadID'			: ";".join(reads),
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
				'IntStart' 	  : row['Position 1'] if row['Chromosome 2'] == 'chrVirus' else row['Position 2'],
				'IntStop' 	  : row['Position 1'] if row['Chromosome 2'] == 'chrVirus' else row['Position 2'],
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
				'IntStart' 	  : row['left_pos'] if 'chr' in row['@left_chr'] else row['right_pos'],
				'IntStop' 	  : row['left_pos'] if 'chr' in row['@left_chr'] else row['right_pos'],
				'Orientation' : 'unknown',
				'ReadID'      : '',
				'id'		  : row_num
				}
			data.append(result)
			row_num += 1
			
	return data	


if __name__ == "__main__":
	main(argv[1:])
