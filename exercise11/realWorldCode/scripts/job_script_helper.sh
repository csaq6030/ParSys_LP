for j in 1 2 4 8
do

	export OMP_NUM_THREADS=$j
	for k in $(seq 10)
	do
		./real
	done

done


