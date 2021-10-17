import sys
import sys, getopt

def main(argv):
	infile = ''

	opts, args = getopt.getopt(argv,"hi:",["ifile="])
	for opt, arg in opts:
		if opt == '-h':
			print('2check_graph.py -i <inputfile>')
			sys.exit()
		elif opt in ("-i", "--ifile"):
			infile = arg

	f = open(infile, 'r')

	line = f.readline()
	n,e = map(int,line.split())

	line_num = 1
	# check each line has exact two number
	while True:
		line = f.readline()
		if not line:
			break
		line_num += 1
		a = line.split()
		if len(a) != 2:
			print(line_num)
			print(a)
			sys.exit()
	f.close()

	print(">> Graph is OK")

if __name__ == "__main__":
   main(sys.argv[1:])
