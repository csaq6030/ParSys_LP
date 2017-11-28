for j in 1 2 4 8
do

	for k in $(seq 7)
	do
		./stencil.o $j 1 100000000 1 2 3
	done

done

