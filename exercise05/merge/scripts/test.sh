#!/bin/sh

for j in  100000 1000000 10000000
do
    for i in $(seq 10)
    do
        $1 $j $2
        echo "Done iteration $i for $j"
    done
done
