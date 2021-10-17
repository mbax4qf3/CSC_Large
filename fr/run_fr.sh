#! /bin/bash -e

echo "*** Graph Preprocessing ***"
echo ""
# unzip dataset wikipedia_link_fr
echo ">> Download fr dataset..."
wget http://konect.cc/files/download.tsv.wikipedia_link_fr.tar.bz2
echo ">> Unzip and clean dataset (wait for minutes)"
tar -jxvf download.tsv.wikipedia_link_fr.tar.bz2
echo "> Copy graph to graph folder..."
cp ./wikipedia_link_fr/out.wikipedia_link_fr ./graph/
echo "> Remove info line from graph..."
sed -i '1d' ./graph/out.wikipedia_link_fr
echo "> Clean graph..."
python3 1graph_clean.py -i ./graph/out.wikipedia_link_fr -o ./graph/1.txt
echo "> Check graph..."
python3 2check_graph.py -i ./graph/1.txt
echo "> Remove tmp files..."
rm ./graph/out.wikipedia_link_fr
echo ">> Dataset extracted and move to graph folder, named 1.txt."

echo ""
echo "*** Build Index ***"
echo ""
# compile
make
echo "compile completed, u_index and u_query are generated."

# build index for original and query
echo "data 1 ori index constructing"
./u_index -g graph/1.txt -l label/1.txt -b n -o degree
echo "data 1 ori querying"
./u_query -l label/1.txt -q query/ -a answer/ -t o
echo "-----------------------------------------------------"

echo "All done."

