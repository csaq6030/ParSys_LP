for j in 8 10 15
do
    for k in $(seq 5)
    do
    mpirun -np $NSLOTS ./a.out $j
    done
done

