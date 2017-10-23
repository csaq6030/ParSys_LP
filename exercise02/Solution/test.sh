#!/bin/sh
#For GNU Time under the gtime command exchange if needed
for i in $(seq 10)
do
    for j in 10 100 500 1000 1500
    do
        #wall time in [s]; user mode cpu; kernel mode; major page fault; minor page fault; max resident set size;avg resident set size
        gtime -f "%e;%U;%S;%F;%R;%M;%t" $1 $j 2>> $1$j.csv
        echo "Done iteration $i for $j"
    done
done
