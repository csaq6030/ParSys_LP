for j in 1 2 4 8
do

	for k in $(seq 3)
	do
		./stencil.o $j 2 512 20 0.5 -0.5 1 0
	done

done

