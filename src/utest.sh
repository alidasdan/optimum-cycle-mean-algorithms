#!/bin/bash

# testing each executable in a tiny sample file; for results on other
# tests, check out all-min-runs.txt and all-max-runs.txt in under
# github/alidasdan/graph-benchmarks.

for i in *.x; 
do 
   ./$i ../input/sample.d -v 1 | grep lambda | awk -v p=$i -v v=1 -v t=40 -v e=0.01 -f utest.awk 
   ./$i ../input/sample.d -v 0 | grep lambda | awk -v p=$i -v v=0 -v t=50 -v e=0.01 -f utest.awk 
done

# EOF


