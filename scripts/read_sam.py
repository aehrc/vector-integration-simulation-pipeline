import pysam
import pdb
from sys import argv

correct = {'chr2-616' : '2', 'chr2-614' : '2', 'chr2-604': '2', 'chr2-576': '2'}
incorrect = {'chr2-612': '2', 'chr2-610': '2', 'chr2-592': '2', 'chr2-580': '2', 'chr2-504' : '1', 'chr2-432' : '1'}


def main(argv):
	samfile = pysam.AlignmentFile(argv[1])
	
	count = 0
	for line in samfile:
		count += 1
		if line.qname in incorrect.keys():
			if line.is_read1 is (incorrect[line.qname] == '1'):
				if line.is_reverse is True:
					ori = 'reverse'
				else:
					ori = 'forward'
				print(f"{line.qname} is {ori}, {line.cigarstring}")
				#print(line.get_aligned_pairs())
				#print()
	
	samfile.close()
	print(f"read {count} lines from samfile {argv[1]}")
	
	
if __name__ == "__main__":
	main(argv)