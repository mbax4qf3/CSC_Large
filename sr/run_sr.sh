#! /bin/bash -e

echo "*** Graph Preprocessing ***"
echo ""
# unzip dataset wikipedia_link_sr
echo ">> Download sr dataset..."
wget http://konect.cc/files/download.tsv.wikipedia_link_sr.tar.bz2
echo ">> Unzip and clean dataset (wait for minutes)"
tar -jxvf download.tsv.wikipedia_link_sr.tar.bz2
echo "> Copy graph to graph folder..."
cp ./wikipedia_link_sr/out.wikipedia_link_sr ./graph/
echo "> Remove info line from graph..."
sed -i '1d' ./graph/out.wikipedia_link_sr
echo "> Clean graph..."
python3 1graph_clean.py -i ./graph/out.wikipedia_link_sr -o ./graph/2.txt
echo "> Check graph..."
python3 2check_graph.py -i ./graph/2.txt
echo "> Remove tmp files..."
rm ./graph/out.wikipedia_link_sr
echo ">> Dataset extracted and move to graph folder, named 2.txt."

echo ""
echo "*** Build Index ***"
echo ""
# compile
make
echo "compile completed, u_index and u_query are generated."

# build index for original and query
echo "data 2 ori index constructing"
./u_index -g graph/2.txt -l label/2.txt -b n -o degree
echo "data 2 ori querying"
./u_query -l label/2.txt -q query/ -a answer/ -t o
echo "-----------------------------------------------------"

echo "All done."

