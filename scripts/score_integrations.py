#!/usr/bin/env python3

# Score simulated vs detected integrations on a per-integration basis

# do this by counting the number of lines in found-info file that 
# fall within a specified window around each integration in sim-info file

from sys import argv
import argparse
import csv
import pdb

# currently supported tools
supported_tools = ('pipeline', 'vifi', 'polyidus')

def main(argv):
	#get arguments
	parser = argparse.ArgumentParser(description='simulate viral integrations')
	parser.add_argument('--sim-info', help='information from simulation', required = True, type=str)
	parser.add_argument('--found-info', help='information from found of simulated reads', required = True, type=str)
	parser.add_argument('--window', help='look this amount either side of simulated hPos for integrations in found', required=False, type=int, default=100)
	parser.add_argument('--output', help='output file', required=False, default='results.tsv')
	parser.add_argument('--summary', help='output summary file', required=False, default='summary.tsv')	
	parser.add_argument('--analysis-tool', choices=supported_tools, required=True)
	
	args = parser.parse_args()
	
	# import information about which integrations were found
	if args.analysis_tool == "pipeline":
		found = parse_results_tsv(args.found_info)
	elif args.analysis_tool == "polyidus":
		found = parse_polyidus(args.found_info)
	elif args.analysis_tool == "vifi":
		found = parse_vifi(args.found_info)
		
	# import infomration about simulated integrations
	sim = parse_tsv(args.sim_info)
	
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
						'hv_count' : 0,
						'vh_count' : 0,
						'total_count': 0,
						'hv_reads' : [],
						'vh_reads' : [],
						'total_reads': []
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
			if overlap(found_start, found_stop, sim_start, sim_stop):
				
				# increment counters
				if result['Orientation'] != "unknown":
					side = result['Orientation']
					sim_results[f"{side}_count"] += 1
					sim_results[f"{side}_reads"].append(result['ReadID'])
				sim_results['total_count'] += 1
				sim_results['total_reads'].append(result['ReadID'])
				
				# mark this integration in found
				if 'found' not in result:
					result['found'] = 1
				else:
					result['found'] += 1
		
		# collapse lists
		for field in ('hv_reads', 'vh_reads', 'total_reads'):
			sim_results[field] = ";".join(sim_results[field])
		
		scored_sim_results.append(sim_results)
	
	scores = {'sim_info'      : args.sim_info,
			  'found_info' : args.found_info,
			  'analysis_tool' : args.analysis_tool,
			  'window'		  : args.window,
			  'tp':0, 'fp':0, 'tn':0, 'fn':0 }
	# check scored_sim_results (which has one entry for each simulated integration) 
	# for true positives and false negatives
	for result in scored_sim_results:
		# integrations that we saw, and expected to see
		if result['total_count'] > 0:
			scores['tp'] += 1
			result['score'] = 'tp'
		# integrations that we didn't
		else:
			scores['fn'] += 1
			result['score'] = 'fn'
	# check found for integrations that we saw but didn't expect to
	fps = {}
	for result in found:
	
		# if this result wasn't in the sim integrations
		if 'found' not in result:
		
			# if this is the first read we've seen from this integration
			id = result['id']
			if result['id'] not in fps:
				scores['fp'] += 1
				
				# we want to output this information - collect these reads together
				fps[id] = {
					'id': id,
					'score': 'fp',
					'chr': result['Chr'],
					'pos': f"{result['IntStart']}/{result['IntStop']}",
					'hv_count' : 0,
					'vh_count' : 0,
					'total_count': 1,
					'hv_reads': [],
					'vh_reads': [],
					'total_reads': [result['ReadID']]
				}
			
			# if it's not the first time
			else:
				fps[id]['total_count'] += 1
				fps[id]['total_reads'].append(result['ReadID'])	
			
			# add orientation-specific info, if available
			if result['Orientation'] == 'hv':
				fps[id]['hv_count'] += 1
				fps[id]['hv_reads'].append(result['ReadID'])
			elif result['Orientation'] == 'vh':
				fps[id]['vh_count'] += 1
				fps[id]['vh_reads'].append(result['ReadID'])		
				
	# add false positives to scored_sim_results for output
	for result_dict in fps.values():
		result_dict['hv_reads'] = ";".join(result_dict['hv_reads'])
		result_dict['vh_reads'] = ";".join(result_dict['vh_reads'])
		result_dict['total_reads'] = ";".join(result_dict['total_reads'])
		scored_sim_results.append(result_dict)
			
	# write results to file
	with open(args.output, "w", newline = '') as outfile:
		# create DictWriter object for scored_sim_results
		output_writer = csv.DictWriter(outfile, delimiter = '\t', 
										fieldnames = ('id', 'score', 'chr', 'pos',
										'hv_count', 'vh_count', 'total_count', 
										'hv_reads', 'vh_reads', 'total_reads'))
		output_writer.writeheader()
		for row in scored_sim_results:
			output_writer.writerow(row)
			

	
	# write summary
	with open(args.summary, "w", newline = '') as outfile:
		# create dictwrite object for summary
		output_writer = csv.DictWriter(outfile, delimiter = '\t',
										fieldnames = ('sim_info', 'found_info',
													   'analysis_tool', 'window',
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
					'Orientation'	: row['Orientation'],
					'ReadID'     	: row['ReadID'],
					'id'			: f"{row['Chr']}/{row['IntStart']}/{row['VirusRef']}/{row['VirusStart']}/{row['ReadID']}"
					}

			data.append(info)
			
	return data
	
def parse_polyidus(path):
	"""
	parse polyidus output file (exactHpvIntegrations.tsv), and produce a data structure 
	similar to parse_results_tsv
	
	assign unique id to each integration, which consists of Chrom/IntegrationSite/ChromVirus/ViralIntegrationSite
	"""
	with open(path, newline = '') as found:
		found_reader = csv.DictReader(found, delimiter = '\t')
		data = []
		
		# get info from rows
		for row in found_reader:
			
			framgents = row['FragmentName'].split(", ")
			
			for fragment in framgents:
				result = {
					'Chr' 		  : row['Chrom'],
					'IntStart' 	  : row['IntegrationSite'],
					'IntStop' 	  : row['IntegrationSite'],
					'Orientation' : 'unknown',
					'ReadID'      : fragment,
					'id'		  : f"{row['Chrom']}/{row['IntegrationSite']}/{row['ChromVirus']}/{row['ViralIntegrationSite']}"
				}
				data.append(result)
			
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
		
		for line in handle:
			# ignore divider lines
			if line[:5] == "##===":
				continue
			# if this is a new integration, clear buffer
			elif line[0] != "#":
				n_reads = int(line.split()[3])
				n_seen = 0
				id = "/".join(line.split())
	
			else:
				# right now we only care about host reads, which always come first

				if n_seen >= n_reads:
					continue
				
				# get info from this line
				info = line.split()
				# strip leading '#' from read ID
				info[0] = info[0][2:]

				ori = 'hv' if info[3] == "True" else 'vh'
				result = {
					'Chr' 		  : info[1],
					'IntStart' 	  : int(info[2]),
					'IntStop' 	  : int(info[2]),
					'Orientation' : ori,
					'ReadID'      : info[0],
					'id'		  : id
				}
				data.append(result)
				n_seen += 1
				
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
	
if __name__ == "__main__":
	main(argv[1:])
