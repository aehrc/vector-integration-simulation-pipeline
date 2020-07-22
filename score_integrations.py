# Score simulated vs detected integrations

from sys import argv
import argparse
import csv
import pdb

def main(argv):
	#get arguments
	parser = argparse.ArgumentParser(description='simulate viral integrations')
	parser.add_argument('--sim-info', help='information from simulation', required = True, type=str)
	parser.add_argument('--analysis-info', help='information from analysis of simulated reads', required = True, type=str)
	parser.add_argument('--window', help='look this amount either side of simulated hPos for integrations in analysis', required=False, type=int, default=1)
	parser.add_argument('--output', help='output file', required=False, default='results.tsv')
	args = parser.parse_args()

	results = []
	
	# import simulated and pipeline information	
	with open(args.sim_info, newline = '') as sim, open(args.analysis_info, newline='') as analysis, open(args.output, "w", newline = '') as outfile:
		
		# create DictReader objects for inputs
		sim_reader = csv.DictReader(sim, delimiter = '\t')
		analysis_reader = csv.DictReader(analysis, delimiter = '\t')
		
		# create DictWriter object for output
		output_writer = csv.DictWriter(outfile, delimiter = '\t', 
										fieldnames = ['id', 'hv_count', 'vh_count', 'hv_reads', 'vh_reads'])
		output_writer.writeheader()
		
		for sim_row in sim_reader:
			# to store the results that we find for each simulated integration
			sim_results = {'id' : sim_row['id'],
							'hv_count' : 0,
							'vh_count' : 0,
							'hv_reads' : [],
							'vh_reads' : []
							}
				
			# look through rows of analysis for matches
			analysis.seek(0)
			for analysis_row in analysis_reader:
				# check for same chromosome in host
				if analysis_row['Chr'] != sim_row['chr']:
					continue
				
				# get start and stop positions to check for overlap in host
				analysis_start = int(analysis_row['IntStart'])
				analysis_stop = int(analysis_row['IntStop'])
				sim_start = int(sim_row['hPos']) - args.window
				sim_stop = int(sim_row['hPos']) + args.window
				if overlap(analysis_start, analysis_stop, sim_start, sim_stop):
					
					# increment counters
					side = analysis_row['Orientation']
					sim_results[f"{side}_count"] += 1
					sim_results[f"{side}_reads"].append(analysis_row['ReadID'])
					
			# collapse lists
			sim_results['hv_reads'] = ";".join(sim_results['hv_reads'])
			sim_results['vh_reads'] = ";".join(sim_results['vh_reads'])
			
			# write row of output
			output_writer.writerow(sim_results)
			
			
	
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