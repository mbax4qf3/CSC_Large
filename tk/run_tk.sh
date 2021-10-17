#! /bin/bash -e

echo "*** Graph Preprocessing ***"
echo ""
# unzip dataset wikipedia_link_tk
echo ">> Download tk dataset..."
wget http://konect.cc/files/download.tsv.wikipedia_link_tk.tar.bz2
echo ">> Unzip and clean dataset (wait for minutes)"
tar -jxvf download.tsv.wikipedia_link_tk.tar.bz2
echo "> Copy graph to graph folder..."
cp ./wikipedia_link_tk/out.wikipedia_link_tk ./graph/
echo "> Remove info line from graph..."
sed -i '1d' ./graph/out.wikipedia_link_tk
echo "> Clean graph..."
python3 1graph_clean.py -i ./graph/out.wikipedia_link_tk -o ./graph/0.txt
echo "> Check graph..."
python3 2check_graph.py -i ./graph/0.txt
echo "> Remove tmp files..."
rm ./graph/out.wikipedia_link_tk
echo ">> Dataset extracted and move to graph folder, named 0.txt."

echo ""
echo "*** Build Index ***"
echo ""
# compile
make
echo "compile completed, u_index and u_query are generated."

# build index for original and query
echo "data 0 ori index constructing"
./u_index -g graph/0.txt -l label/0.txt -b n -o degree
echo "data 0 ori querying"
./u_query -l label/0.txt -q query/ -a answer/ -t o
echo "-----------------------------------------------------"

echo "All done."

