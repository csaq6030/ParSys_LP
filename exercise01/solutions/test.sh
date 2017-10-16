#!/bin/sh
#For GNU Time under the gtime command exchange if needed
for i in $(seq 10)
do
    for j in 10 100 500 1000 1500
    do
        gtime -f "%e;%U;%S" ./mmul $j 2>> mmul_$j.csv
        echo "Done iteration $i for $j"
    done
done
