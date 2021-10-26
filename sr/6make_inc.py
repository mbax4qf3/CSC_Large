from collections import defaultdict
import random

inum = '2'
infile = "./graph/" + inum + ".txt"
outfilei = "./graph/" + inum + "_inc.txt"

f = open(infile, 'r')
f_outi = open(outfilei, 'w')
f_outi.write("300\n")

line = f.readline()
node_num, edge_num = map(int,line.split())
nodes = set()
e = defaultdict(int)

while True:
	line = f.readline()
	if not line:
		break
	a,b = map(int,line.split())
	if a not in nodes:
		nodes.add(a)
	if b not in nodes:
		nodes.add(b)
	e[(a,b)] = 1

inc_cout = 0
while (inc_cout < 300):
	a = random.randint(1,node_num-1)
	b = random.randint(1,node_num-1)
	if a not in nodes or b not in nodes:
		continue
	if (e[(a,b)]==1 or e[(b,a)==1]):
		continue
	e[(a,b)] = 1
	str_out = str(a) + " " + str(b+node_num) + "\n"
	f_outi.write(str_out)
	inc_cout+=1

