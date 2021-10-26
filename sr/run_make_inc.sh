#! /bin/bash -e

echo "Make bipartite and inc file"
python3 6make_inc.py

make
echo "compiled completed"
echo "data 2_bi edge adding and querying"
./u_query -l label/2_bi.txt -i graph/2_inc.txt -q query/ -a answer/ -t i

echo "Done"
