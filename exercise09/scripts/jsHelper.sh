for k in $(seq 5)
do
	mpirun -np $NSLOTS ./stencil512.o 
done

for k in $(seq 5)
do
	mpirun -np $NSLOTS ./stencil768.o 
done

