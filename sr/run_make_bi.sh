#! /bin/bash -e

echo "Make bipartite and inc file..."
python3 5make_bi_inc.py
sed -i '1s/^/6350020 142760498\n/' ./graph/2_bi.txt
echo "Check file..."
python3 2check_graph.py -i ./graph/2.txt

make
echo "compiled completed"

echo "data 2_bi index constructing"
./u_index -g graph/2_bi.txt -l label/2_bi.txt -b y -o bidegree
echo "data 2_bi edge adding and querying"
./u_query -l label/2_bi.txt -i graph/2_inc.txt -q query/ -a answer/ -t i

echo "Done"
