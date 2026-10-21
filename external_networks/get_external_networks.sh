#!/bin/sh
# Downloads the ten public edge lists used in Section 7.5 into ./raw (Stanford Network Analysis Project, Leskovec & Krevl 2014).
# The files are treated as untrusted data: the analysis reads them only as two integer columns.
mkdir -p raw && cd raw || exit 1
B=https://snap.stanford.edu/data
for f in ca-GrQc ca-HepTh ca-CondMat ca-AstroPh email-Enron facebook_combined loc-brightkite_edges loc-gowalla_edges; do curl -fsSL -o "$f.txt.gz" "$B/$f.txt.gz"; done
curl -fsSL -o com-amazon.txt.gz "$B/bigdata/communities/com-amazon.ungraph.txt.gz"
curl -fsSL -o com-dblp.txt.gz   "$B/bigdata/communities/com-dblp.ungraph.txt.gz"
ls -l
