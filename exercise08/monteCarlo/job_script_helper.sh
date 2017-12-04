for i in 1000000 10000000 100000000
do

	for j in 1 2 4 8
	do
	
		export OMP_NUM_THREADS=$j
		for k in $(seq 10)
		do
			./monteCarlo.o $i par
		done

	done

done
