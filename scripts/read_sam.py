import pysam
import pdb
from sys import argv

correct = {'chr2-612' : '2', 'chr2-610' : '2', 'chr2-604' : '2', 'chr2-582' : '1', 'chr2-580' : '2'}
incorrect = {'chr2-616' : '2', 'chr2-614' : '2', 'chr2-602' : '2', 'chr2-592' : '2', 'chr2-586' : '1'}


def main(argv):
	samfile = pysam.AlignmentFile(argv[1])
	
	count = 0
	for line in samfile:
		count += 1
		if line.qname in correct.keys():
			if line.is_reverse is True:
				ori = 'reverse'
			else:
				ori = 'forward'
			if line.is_read1 is True:
				read_num = '1'
			else:
				read_num = '2'
				
			if read_num != correct[line.qname]:
				continue

			print(f"{line.qname}/{read_num} is {ori} with cigar {line.cigarstring} and pos {line.reference_start}")
				#print(line.get_aligned_pairs())
				#print()
	
	samfile.close()
	print(f"read {count} lines from samfile {argv[1]}")
	
	
if __name__ == "__main__":
	main(argv)