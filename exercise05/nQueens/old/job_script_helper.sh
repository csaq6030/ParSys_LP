for i in 8 10 15
do

	for j in 4 8
	do
	
		export OMP_NUM_THREADS=$j
		for k in $(seq 10)
		do
			./nQueens.o $i 
		done

	done

done
