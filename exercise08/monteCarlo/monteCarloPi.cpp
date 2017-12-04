#include <iostream>
#include <random>
#include <cmath>
#include <omp.h>
#include <mpi.h>

#include "chrono_timer.h"

using namespace std;

double MCS_Seq(unsigned long n){
	unsigned long c = 0;
	double x, y;
	random_device rd;	// random seed
	mt19937 gen(rd()); 	// random generator
	uniform_real_distribution<double> dis(0, 1.0);
	for(unsigned i = 0; i < n; i++){
		x = dis(gen);
		y = dis(gen);
		if(sqrt(x*x + y*y) <= 1)
			c++;
	}
	return (double) c/n*4.0;
}



int main(int argc, char** argv) {

	if(argc != 2) return EXIT_FAILURE;
	unsigned long long n = atoi(argv[1]);
	if(n == 0) return EXIT_FAILURE;

	unsigned long c = 0;
	double pi;
	int myid;
	int size;
	int subIter;
	int reducedCount;
	double x, y;
	// start(std::chrono::high_resolution_clock::now());
	MPI_Init(NULL, NULL);				
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);	
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	subIter = n/(unsigned long)size;

	random_device rd;
	mt19937 gen(rd() + myid); 
	uniform_real_distribution<double> dis(0, 1.0);
	int i;
	if(myid == 0){
		for (i=0; i<subIter + n%size; ++i){
			x = dis(gen);
			y = dis(gen);
			if(sqrt(x*x + y*y) <= 1)
				c++;	
		}	
	}else{
		for (i=0; i<subIter; ++i){
			x = dis(gen);
			y = dis(gen);
			if(sqrt(x*x + y*y) <= 1)
				c++;	
		}
	}
	cout << "procces id: " << myid << "   iter: " << i << endl;
	MPI_Reduce(&c, &reducedCount, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Finalize(); 
	
	if (myid == 0){
		pi = ((double)reducedCount/(double)n)*4.0;		
		cout.precision(11);
		cout << fixed;
		cout<< "Estimated value of pi is: " << pi << endl;
		cout<< "Actual    Value of pi is: 3.14159265359" << endl;

		//auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start);
		//cout << elapsed.count() << std::endl;
	}


	return 0;
}









