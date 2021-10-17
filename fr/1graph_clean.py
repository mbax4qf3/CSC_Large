from collections import defaultdict
import sys, getopt

def main(argv):
	infile = ''
	outfile = ''

	opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
	for opt, arg in opts:
		if opt == '-h':
			print('1graph_clean.py -i <inputfile> -o <outputfile>')
			sys.exit()
		elif opt in ("-i", "--ifile"):
			infile = arg
		elif opt in ("-o", "--ofile"):
			outfile = arg
	f = open(infile, 'r')
	f_out = open(outfile, 'w')

	maxNode = 0
	num_edge = 0

	split_char = " "
	e = defaultdict(int)
	while True:
		line = f.readline().strip()
		if not line:
			break
	
		a,b = line.split()
		a = int(a)
		b = int(b)
		# check self loop
		if (a == b):
			continue

		# check multiple edge
		if (e[(a,b)] == 1):
			continue
		e[(a,b)] = 1
		num_edge += 1

		maxNode = max(maxNode, a, b)

		str_out = str(a) + split_char + str(b) + "\n"
		f_out.write(str_out)
	f.close()
	f_out.close()

	# write N, M info
	head = str(maxNode+1) + split_char + str(num_edge) + "\n"
	with open(outfile, 'r+') as f_out2:
		content = f_out2.read()
		f_out2.seek(0, 0)
		f_out2.write(head + content)

if __name__ == "__main__":
   main(sys.argv[1:])
