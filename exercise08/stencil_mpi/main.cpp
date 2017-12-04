#include <iostream>
#include <cstdlib>
#include <chrono>
#include <mpi.h>
#include <cmath>
#include <algorithm>

#include "wrapper.h"

using namespace std;

void printError() {
    cout << "Usage: [number of threads] [dimension] [size] [epsilon] [left] [right] [up] [down] [front] [back]" << endl
    << "up and down are needed only for 2D and 3D" << endl
    << "front and back only for 3D" << endl
    << "num of number of threads is ignored, just input something in this place" << endl;
}


int main(int argc, char* argv[]) {
    int dim = 0;
    int size = 0;
    float epsilon = 0.;
    float left = 0;
    float right = 0;
	float reducedProgress = 0;
    
    int myid, worldSize;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);

	
    if ( argc > 11  || argc < 7) {
        printError();
		MPI_Finalize(); 
        return EXIT_FAILURE;
    }
    

    
    if (argc >= 7) {
		//omp_set_num_threads(atoi(argv[1]));

        dim = atoi(argv[2]);
        size = atoi(argv[3]);
        epsilon = atof(argv[4]);
        left = atof(argv[5]);
        right = atof(argv[6]);
    }

	const int blockSize = size/(int)worldSize;
	
    //neighbour indexing
	int leftid, rightid;
	if(myid!=0){
		leftid = myid -1;
	}
    
	if(myid!=worldSize-1){
		rightid = myid +1;
	}
    
    if (atoi(argv[2]) == 1 && argc == 7) { //1D

		if(myid == 0){
			Wrapper firstBlock(1, blockSize + 1);
			//boundary left
			firstBlock.arrayA[0] = left;
			firstBlock.arrayB[0] = left;
			//cout << "id: " << myid << endl;
			while(true){
				float progress = 0;
				for (int i=1; i<blockSize + 2; ++i){
					firstBlock.arrayB[i] = (firstBlock.arrayA[i-1] + firstBlock.arrayA[i] + firstBlock.arrayA[i + 1])/3; 
					//progress += std::abs(myBlock.arrayB[i]-myBlock.arrayA[i]);
				}
				for (int i=1; i<blockSize + 1; ++i){
					firstBlock.arrayA[i] = (firstBlock.arrayB[i-1] + firstBlock.arrayB[i] + firstBlock.arrayB[i + 1])/3;
					progress += std::abs(firstBlock.arrayB[i]-firstBlock.arrayA[i]);
				}
				MPI_Allreduce(&progress, &reducedProgress, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
				if(reducedProgress < epsilon){
					//send result - none master //collect data
					float *result = new float[size];
					for(int i = 1; i< blockSize+1; i++){
						result[i-1] = firstBlock.arrayA[i];
					}
					for(int i = 1; i< worldSize-1; i++){
						MPI_Recv(&result[i*blockSize], blockSize, MPI_FLOAT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					}
					MPI_Recv(&result[(worldSize-1)*blockSize], blockSize+size%worldSize, MPI_FLOAT, worldSize-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					for(int i = 0; i< size; i++){
						cout << result[i] << "  ";
					}
					cout << endl;
					delete[] result;
					break;
				}
				MPI_Send(&firstBlock.arrayA[blockSize], 2, MPI_FLOAT, rightid, 0, MPI_COMM_WORLD);
				MPI_Recv(&firstBlock.arrayA[blockSize + 1], 2, MPI_FLOAT, rightid, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}	
			//
		}else if(myid == (worldSize-1) ){
			Wrapper lastBlock(1, blockSize + size % worldSize + 1);
			//boundary right
			lastBlock.arrayA[blockSize + size % worldSize+ 2] = right;
			lastBlock.arrayB[blockSize + size % worldSize+ 2] = right;
			
				//cout << "id: " << myid << endl;
			while(true){
				float progress = 0;
				for (int i=1; i<blockSize + size % worldSize + 2; ++i){
					lastBlock.arrayB[i] = (lastBlock.arrayA[i - 1] + lastBlock.arrayA[i] + lastBlock.arrayA[i + 1]) / 3;
					//progress += std::abs(myBlock.arrayB[i]-myBlock.arrayA[i]);
				}
				for (int i=2; i<blockSize + size % worldSize + 2; ++i){
					lastBlock.arrayA[i] = (lastBlock.arrayB[i - 1] + lastBlock.arrayB[i] + lastBlock.arrayB[i + 1]) / 3;
					progress += std::abs(lastBlock.arrayB[i] - lastBlock.arrayA[i]);
				}
				MPI_Allreduce(&progress, &reducedProgress, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
				if(reducedProgress < epsilon){
					//send result
					cout << "blocksize: " << blockSize << endl;
					for(int i = 0; i< blockSize + size % worldSize + 3; i++)
						cout << lastBlock.arrayA[i] << " ";
					cout << endl;
					MPI_Send(&lastBlock.arrayA[2], blockSize + size % worldSize, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
					break;
				}
				//cout << "id: " << myid << endl;
				MPI_Send(&lastBlock.arrayA[2], 2, MPI_FLOAT, leftid, 0, MPI_COMM_WORLD);
				MPI_Recv(&lastBlock.arrayA[0], 2, MPI_FLOAT, leftid, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				
			}	
		}else {
			Wrapper myBlock(1, blockSize + 2);
			
			while(true){
				float progress = 0;
				for (int i = 1; i< blockSize - 1; ++i){
					myBlock.arrayB[i] = (myBlock.arrayA[i - 1] + myBlock.arrayA[i] + myBlock.arrayA[i + 1]) / 3;
					//progress += std::abs(myBlock.arrayB[i]-myBlock.arrayA[i]);
				}
				for (int i=2; i < blockSize - 2; ++i){
					myBlock.arrayA[i] = (myBlock.arrayB[i - 1] + myBlock.arrayB[i] + myBlock.arrayB[i + 1]) / 3;
					progress += std::abs(myBlock.arrayB[i]-myBlock.arrayA[i]);
				}
				MPI_Allreduce(&progress, &reducedProgress, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
				if(reducedProgress < epsilon){
					//cout << "id: " << myid << endl;
					//send result
					MPI_Send(&myBlock.arrayA[2], blockSize, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
					break;
				}
				MPI_Send(&myBlock.arrayA[2], 2, MPI_FLOAT, leftid, 0, MPI_COMM_WORLD);
				MPI_Send(&myBlock.arrayA[blockSize], 2, MPI_FLOAT, rightid, 0, MPI_COMM_WORLD);
				MPI_Recv(&myBlock.arrayA[0], 2, MPI_FLOAT, leftid, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Recv(&myBlock.arrayA[blockSize + 2], 2, MPI_FLOAT, rightid, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
		}
		
        
    } else if (atoi(argv[2]) == 2 && argc == 9) { //2D
        const float up = atof(argv[7]);
        const float down = atof(argv[8]);

        
    } else if (atoi(argv[2]) == 3 && argc == 11) { //3D
        const float up = atof(argv[7]);
        const float down = atof(argv[8]);
        const float front = atof(argv[9]);
        const float back = atof(argv[10]);

        
        //cout <<"num_t: " << num_threads << " dim: " << dim << " size: " << size << " epsilon: " << epsilon << " left: " << left << " right: " << right << " up: " << up << " down: " << down << " front: " << front << " back: " << back << endl;
        
    } else {
        printError();
		MPI_Finalize(); 
        return EXIT_FAILURE;
    }
	cout << "destroy: " << myid << endl;
	MPI_Finalize();
    
    return EXIT_SUCCESS;
}
