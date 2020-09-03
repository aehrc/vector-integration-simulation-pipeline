import pysam
import pdb
from sys import argv

correct = {'chr2-594', 'chr2-578'}
incorrect = {'chr2-574', 'chr2-578', 'chr2-530', 'chr2-376', 'chr2-268'}


def main(argv):
	samfile = pysam.AlignmentFile(argv[1])
	
	count = 0
	for line in samfile:
		count += 1
		if line.qname in incorrect:
			if line.is_reverse is True:
				ori = 'reverse'
			else:
				ori = 'forward'
			if line.is_read1 is True:
				read_num = '1'
			else:
				read_num = '2'

			print(f"{line.qname}/{read_num} is {ori} with cigar {line.cigarstring} and pos {line.reference_start}")
				#print(line.get_aligned_pairs())
				#print()
	
	samfile.close()
	print(f"read {count} lines from samfile {argv[1]}")
	
	
if __name__ == "__main__":
	main(argv)