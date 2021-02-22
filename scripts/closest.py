#!/ur/bin/env python3

# for each interval in an 'a' file, find the closest integration in a 'b' file
# similar to bedtools 'closest', but use different kinds of distances

# types of distances:
## shortest: bedtools-style, the shortest distance between the starts and stops
## coords_mean: mean of absolute differences between starts and stops
## coords_min: min of absolute differences between starts and stops
## midpoint: absolute difference between midpoints

import argparse
import sys
import csv
import pdb

def main(args):
	parser = argparse.ArgumentParser(description='Apply threshold to distances to calculate TP, FP, FN')
	parser.add_argument('-a', help='A file with intervals - f', required=True)
	parser.add_argument('-b', help='B file - find closest interval in this file for each interval in a file', required=True)
	parser.add_argument('-o', help='Output file')	
	args = parser.parse_args(args[1:])
	
	a_bed = read_bed(args.a)
	b_bed = read_bed(args.b)
	
	for i_a in a_bed:
		for dist in dist_types:
			
			i_a += list(find_closest(i_a, b_bed, dist))
			
	header_base = ('closest_chr', 'closest_start', 'closest_stop', 'closest_ori', 'd')
	header = ['chr', 'start', 'stop', 'ori']
	for dist in dist_types:
		header += [f"{dist.__qualname__}_{col}" for col in header_base]
	header = "\t".join(header)
	
	# if no output file, print to stdout
	if args.o is None:
		print(header)
	else:
		outfile = open(args.o, 'w')
		outfile.write(header + "\n")
		
	for line in a_bed:
		write_line = "\t".join([str(i) for i in line])
		if args.o is None:
			print(write_line)
		else:
			outfile.write(write_line + "\n")

	if args.o is not None:
		outfile.close()
	
	
def find_closest(i_a, intervals_b, dist_method):
	d = 1e12
	closest = (".", -1, -1, ".")
	for i_b in intervals_b:
	
		# if chromosomes don't match
		if i_a[0] != i_b[0]:
			continue
			
		# if orientations don't match
		if i_a[3] != i_b[3]:
			continue		
			
		# check distance
		check_coords(i_a[1], i_a[2], i_b[1], i_b[2])
		dist = dist_method(i_a[1], i_a[2], i_b[1], i_b[2])
		if dist < d:
			d = dist
			closest = i_b
			
	if closest == (".", -1, -1, "."):
		d = -1
			
	return *closest, d

def check_coords(start_1, stop_1, start_2, stop_2):
	assert stop_1 >= start_1
	assert stop_2 >= start_2
	
def shortest(start_1, stop_1, start_2, stop_2):
	"""
	find the shortest distance between two coordinates.  If they overlap, the distance is 0
	"""

	# interval 1 competley to the left of interval 2
	if stop_1 < start_2:
		return start_2 - stop_1
		
	# interval 2 completely to the left of interval 1
	if stop_2 < start_1:
		return start_1 - stop_2
		
	# overlap
	return 0
	
def coords_mean(start_1, stop_1, start_2, stop_2):	
	d1 = abs(start_1 - start_2)
	d2 = abs(stop_1 - stop_2)
	return((d1 + d2) / 2)
	
def coords_min(start_1, stop_1, start_2, stop_2):	
	d1 = abs(start_1 - start_2)
	d2 = abs(stop_1 - stop_2)
	return(min(d1, d2))
	
def midpoint(start_1, stop_1, start_2, stop_2):	
	d1 = (stop_1 - start_1)/2 + start_1
	d2 = (stop_2 - start_2)/2 + start_2
	return(abs(d2 - d1))
	
def read_bed(file):
	intervals = []
	id = 0
	with open(file, 'r') as bedfile:
		bed = csv.reader(bedfile, delimiter='\t')
		for line in bed:
			# check for integer start and stop
			try:
				start, stop = int(line[1]), int(line[2])
			except ValueError:
				raise ValueError(f"Input file {file} does not look like a bed file: line {line} contains non-integer coordinates")
			
			# check for start not larger than stop
			if start > stop:
				raise ValueError(f"Input file {file} does not look like a bed file: line {line} start larger than stop")
			
			# check for orientation
			try:
				ori = line[3]
			except IndexError:
				ori = '+'
			
			intervals.append([line[0], start, stop, ori])	
			id += 1
			
	return intervals
	
if __name__ == "__main__":
	dist_types = [shortest, coords_mean, coords_min, midpoint]
	main(sys.argv)
