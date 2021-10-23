#! /bin/bash -e

echo "*** Graph Preprocessing ***"
echo ""
echo ">> Download sr dataset..."
wget http://konect.cc/files/download.tsv.zhishi-hudong-relatedpages.tar.bz2
echo ">> Unzip and clean dataset (wait for minutes)"
tar -jxvf download.tsv.zhishi-hudong-relatedpages.tar.bz2
echo "> Remove info line from graph..."
sed -i '1,2d' ./zhishi-hudong-relatedpages/out.zhishi-hudong-relatedpages
echo "> Clean graph..."
python3 1graph_clean.py -i ./zhishi-hudong-relatedpages/out.zhishi-hudong-relatedpages -o 1.txt
echo "> Check graph..."
python3 2check_graph.py -i 1.txt
echo "> Remove tmp files..."
rm -r ./zhishi-hudong-relatedpages
rm download.tsv.zhishi-hudong-relatedpages.tar.bz2
echo ">> Dataset extracted, named 1.txt."
