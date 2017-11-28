for j in 1 2 4 8
do

	for k in $(seq 7)
	do
		./stencil.o $j 3 40 2 2 2 2 2 2 2
	done

done

