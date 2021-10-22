from collections import defaultdict
import random

inum = '2'
infile = "./graph/" + inum + ".txt"
outfile = "./graph/" + inum + "_bi.txt"
outfilei = "./graph/" + inum + "_inc.txt"

f = open(infile, 'r')
f_out = open(outfile, 'w')
f_outi = open(outfilei, 'w')
f_outi.write("200\n")

line = f.readline()
node_num, edge_num = map(int,line.split())
nodes = set()
e = defaultdict(int)
max_n = 0
while True:
	line = f.readline()
	if not line:
		break
	a,b = map(int,line.split())
	if a not in nodes:
		nodes.add(a)
		str_out = str(a+node_num) + " " + str(a) + "\n"
		f_out.write(str_out)
		edge_num += 1
		max_n = max(max_n,a)
	if b not in nodes:
		nodes.add(b)
		str_out = str(b+node_num) + " " + str(b) + "\n"
		f_out.write(str_out)
		edge_num += 1
		max_n = max(max_n,b)
	e[(a,b)] = 1
	str_out = str(a) + " " + str(b+node_num) + "\n"
	f_out.write(str_out)

print(str(max_n+node_num) + " " + str(edge_num))

inc_cout = 0
while (inc_cout < 200):
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

