
do

	export OMP_NUM_THREADS = 8
	for k in $(seq 1)
	do
		./a.out
	done

done


