for j in 1 2 4 8
do
    echo "run with $j"
	for k in $(seq 10)
	do
		./stencil.o $j 2 512 20 0.5 -0.5 1 0
	done
done

