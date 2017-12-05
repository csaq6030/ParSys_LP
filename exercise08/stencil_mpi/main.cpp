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

	
	
    //neighbour indexing
	int leftid, rightid;
	if(myid!=0){
		leftid = myid -1;
	}
    
	if(myid!=worldSize-1){
		rightid = myid +1;
	}
    
    const int blockSize = size/(int)worldSize;
    
    if (atoi(argv[2]) == 1 && argc == 7) { //1D
        
        
		if(myid == 0){
			Wrapper firstBlock(1, blockSize + 1); // wrapper has internal +2
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
			Wrapper lastBlock(1, blockSize + size % worldSize + 1); // wrapper has internal +2
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
			Wrapper myBlock(1, blockSize + 2); // wrapper has internal +2
			
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
        
        
        if(myid == 0){ // first stripe
            const int n = size + 2; //m equals i in A[i][j]
            const int m = (blockSize + 3); //m equals j in A[i][j]
            
            float *arrayA = new float[n * m]();
            float *arrayB = new float[n * m]();
            //boundary left
            for (int i = 0; i < n; i++) {
                arrayA[i * m] = left;
                arrayB[i * m] = left;
            }
            
            //boundary up and down
            for (int i = 0; i < m; i++) {
                arrayA[i] = up;
                arrayB[i] = up;
                
                arrayA[(n-1) * m + i] = down; //works because m - n == 1 and the last element is never used
                arrayB[(n-1) * m + i] = down;
            }
            
            cout << "done init first" << endl;
            float *buffer = new float[size * 2]();
            
            while(true){
                float progress = 0;
                for (int i= 1; i < n - 1; i++){
                    for (int j = 1; j < m - 1; j++) {
                    arrayB[i * m + j] = (arrayA[i * m + j - 1] + arrayA[i * m + j] + arrayA[i * m + j + 1] +
                                         arrayA[(i - 1) * m + j] + arrayA[(i + 1) * m + j])/5;
                    }
                }
                
                for (int i= 1; i < n - 1; i++){
                    for (int j = 1; j < m - 2; j++) {
                        arrayA[i * m + j] = (arrayB[i * m + j - 1] + arrayB[i * m + j] + arrayB[i * m + j + 1] +
                                             arrayB[(i - 1) * m + j] + arrayB[(i + 1) * m + j])/5;
                        progress += std::abs(arrayB[i * m + j] - arrayA[i * m + j]);
                    }
                }
                //cout << "iter first" << endl;
                MPI_Allreduce(&progress, &reducedProgress, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
                if(reducedProgress < epsilon){
                    //send result - none master //collect data
                    //own
                    float *result = new float[size * size]();
                    for(int i = 1; i < size + 1; i++){
                        for (int j = 1; j < blockSize + 1; j++) {
                            result[(i - 1) * size + (j - 1)] = arrayA[i * m + j];
                        }
                    }
                    //other
                    float *tmp1 = new float[size * blockSize]();
                    for(int b = 1; b < worldSize - 1; b++){
                        MPI_Recv(tmp1, size * blockSize, MPI_FLOAT, b, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        for(int i = 0; i < size; i++){
                            for (int j = 0; j < blockSize; j++) {
                                result[i * size + j + b * blockSize] = tmp1[i * blockSize + j];
                            }
                        }
                    }
                    delete[] tmp1;
                    
                    float *tmp2 = new float[size * blockSize + size % worldSize]();
                    //last result
                    MPI_Recv(tmp2,size * (blockSize + size % worldSize), MPI_FLOAT, worldSize-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    for(int i = 0; i < size ; i++){
                        for (int j = 0; j < blockSize + size % worldSize; j++) {
                            result[i * size + j + (worldSize - 1) * blockSize] = tmp2[i * (blockSize + size % worldSize) + j];
                            //result[i * size + j + (worldSize - 1) * blockSize] = 1111; //little test
                        }
                    }
                    delete[] tmp2;
                    
                    //Print result
                    for(int i = 0; i < size; i++){
                        for (int j = 0; j < size; j++) {
                            cout << result[i * size + j] << "  ";
                        }
                        cout << endl;
                    }
                    cout << endl;
                    
                    delete[] buffer;
                    delete[] result;
                    delete[] arrayA;
                    delete[] arrayB;
                    break;
                }
                
                //fill buffer
                for (int i = 1; i < size + 1; i++) {
                    for (int j = blockSize - 1; j < blockSize + 1; j++) {
                        buffer[(i - 1) * 2 + j - blockSize + 1] = arrayA[i * m + j];
                    }
                }
                
                MPI_Send(buffer, size * 2, MPI_FLOAT, rightid, 0, MPI_COMM_WORLD);
                MPI_Recv(buffer, size * 2, MPI_FLOAT, rightid, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                
                /* //Print buffer
                for (int i = 0; i < size; i++) {
                    for (int j = 0; j < 2; j++) {
                        cout << buffer[i * 2 + j] << " ";
                    }
                    cout << endl;
                }
                cout << endl; */
                
                //read buffer
                for (int i = 1; i < size + 1; i++) {
                    for (int j = 0; j < 2; j++) {
                        arrayA[i * m + j + blockSize + 1] = buffer[(i - 1) * 2 + j];
                    }
                }
                
                
            }
        }else if(myid == (worldSize-1) ){ //last stripe
            const int n = size + 2; //m equals i in A[i][j]
            const int m = (blockSize + 3 + size % worldSize); //m equals j in A[i][j]
            
            float *arrayA = new float[n * m]();
            float *arrayB = new float[n * m]();
            
            //boundary right
            for (int i = 0; i < n; i++) {
                arrayA[i * m + m - 1] = right;
                arrayB[i * m  + m - 1] = right;
            }
            
            //boundary up and down
            for (int i = 0; i < m; i++) {
                arrayA[i] = up;
                arrayB[i] = up;
                
                arrayA[(n-1) * m + i] = down;
                arrayB[(n-1) * m + i] = down;
            }
            
            /*
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < m; j++) {
                    cout << arrayA[i * m + j] << " ";
                }
                cout << endl;
            }
            cout << endl; */
            
            float *buffer = new float[size * 2]();
            
            cout << "done init last" << endl;
            
            while(true){
                float progress = 0;
                for (int i = 1; i < n - 1; i++){
                    for (int j = 1; j < m - 1; j++) {
                        arrayB[i * m + j] = (arrayA[i * m + j - 1] + arrayA[i * m + j] + arrayA[i * m + j + 1] +
                                             arrayA[(i - 1) * m + j] + arrayA[(i + 1) * m + j])/5;
                    }
                }
                
                for (int i = 1; i < n - 1; i++){
                    for (int j = 2; j < m - 1; j++) {
                        arrayA[i * m + j] = (arrayB[i * m + j - 1] + arrayB[i * m + j] + arrayB[i * m + j + 1] +
                                             arrayB[(i - 1) * m + j] + arrayB[(i + 1) * m + j])/5;
                        progress += std::abs(arrayB[i * m + j] - arrayA[i * m + j]);
                    }
                }
                //cout << "iter last" << endl;
                MPI_Allreduce(&progress, &reducedProgress, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
                if(reducedProgress < epsilon){
                    cout << "result last sending" << endl;
                    //send result
                    float *tmp = new float[size * (blockSize + size % worldSize)]();
                    
                    for (int i = 1; i < size + 1; i++) {
                        for (int j = 2; j < m - 1; j++) {
                            tmp[(i - 1) * (blockSize + size % worldSize) + (j - 2)] = arrayA[i * m + j];
                        }
                    }
                    
                    /*
                    for (int i = 0; i < n; i++) {
                        for (int j = 0; j < m; j++) {
                            cout << arrayA[i * m + j] << " ";
                        }
                        cout << endl;
                    }
                    cout << endl;
                    
                    for (int i = 0; i < size; i++) {
                        for (int j = 0; j < blockSize + size % worldSize; j++) {
                            cout << tmp[i * (blockSize + size % worldSize) + j] << " ";
                        }
                        cout << endl;
                    }
                    cout << endl; */
                    
                    
                    MPI_Send(tmp, (blockSize + size % worldSize) * size, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
                    
                    delete[] tmp;
                    delete[] buffer;
                    delete[] arrayA;
                    delete[] arrayB;
                    break;
                }
                //cout << "fill buffer last" << endl;
                //fill buffer
                for (int i = 1; i < size + 1; i++) {
                    for (int j = 2; j < 4; j++) {
                        buffer[(i - 1) * 2 + j - 2] = arrayA[i * m + j];
                    }
                }
                //cout << "done fill buffer last, now sending" << endl;
                MPI_Send(buffer, size * 2, MPI_FLOAT, leftid, 0, MPI_COMM_WORLD);
                //cout << "send last, now recv" << endl;
                MPI_Recv(buffer, size * 2, MPI_FLOAT, leftid, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                //cout << "done recv last, now reading buffer" << endl;
                
                //print buffer
                /*
                for (int i = 0; i < size; i++) {
                    for (int j = 0; j < 2; j++) {
                        cout << buffer[i * 2 + j] << " ";
                    }
                    cout << endl;
                } */
                
                //read buffer
                for (int i = 1; i < size + 1; i++) {
                    for (int j = 0; j < 2; j++) {
                        //cout << "i: " << i << " " << "j: " << j << endl;
                        arrayA[i * m + j] = buffer[(i - 1) * 2 + j];
                    }
                }
                
                //cout << "done reading buffer last" << endl;
                
            }
        }else {
            const int n = size + 2; //m equals i in A[i][j]
            const int m = (blockSize + 4); //m equals j in A[i][j]
            
            float *arrayA = new float[n * m]();
            float *arrayB = new float[n * m]();
            
            //boundary up and down
            for (int i = 0; i < m; i++) {
                arrayA[i] = up;
                arrayB[i] = up;
                
                arrayA[(n-1) * m + i] = down; //works because m - n == 1 and the last element is never used
                arrayB[(n-1) * m + i] = down;
            }
            
            float *bufferL = new float[size * 2]();
            float *bufferR = new float[size * 2]();
            
            while(true){
                float progress = 0;
                for (int i = 1; i < n - 1; i++){
                    for (int j = 1; j < m - 1; j++) {
                        arrayB[i * m + j] = (arrayA[i * m + j - 1] + arrayA[i * m + j] + arrayA[i * m + j + 1] +
                                             arrayA[(i - 1) * m + j] + arrayA[(i + 1) * m + j])/5;
                    }
                }
                
                for (int i = 1; i < n - 1; i++){
                    for (int j = 2; j < m - 2; j++) {
                        arrayA[i * m + j] = (arrayB[i * m + j - 1] + arrayB[i * m + j] + arrayB[i * m + j + 1] +
                                             arrayB[(i - 1) * m + j] + arrayB[(i + 1) * m + j])/5;
                        progress += std::abs(arrayB[i * m + j] - arrayA[i * m + j]);
                    }
                }
                MPI_Allreduce(&progress, &reducedProgress, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
                if(reducedProgress < epsilon){
                    //send result
                    float *tmp = new float[size * blockSize]();
                    
                    for (int i = 1; i < n - 1; i++) {
                        for (int j = 2; j < m - 1; j++) {
                            tmp[(i - 1) * blockSize + (j - 2)] = arrayA[i * m + j];
                        }
                    }
                    
                    for (int i = 0; i < n; i++) {
                        for (int j = 0; j < m; j++) {
                            cout << arrayA[i * m + j] << " ";
                        }
                        cout << endl;
                    }
                    cout << endl;
                    
                    for (int i = 0; i < size; i++) {
                        for (int j = 0; j < blockSize ; j++) {
                            cout << tmp[i * blockSize + j] << " ";
                        }
                        cout << endl;
                    }
                    cout << endl;
                    
                    MPI_Send(tmp, blockSize * size, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
                    
                    delete[] tmp;
                    delete[] bufferL;
                    delete[] bufferR;
                    delete[] arrayA;
                    delete[] arrayB;
                    break;
                }
                
                //fill buffer
                for (int i = 1; i < size + 1; i++) {
                    for (int j = 2; j < 4; j++) {
                        bufferL[(i - 1) * 2 + j - 2] = arrayA[i * m + j];
                    }
                }
                
                for (int i = 1; i < size + 1; i++) {
                    for (int j = blockSize; j < blockSize + 2; j++) {
                        bufferR[(i - 1) * 2 + j - blockSize] = arrayA[i * m + j];
                    }
                }
                
                MPI_Send(bufferL, size * 2, MPI_FLOAT, leftid, 0, MPI_COMM_WORLD);
                MPI_Send(bufferR, size * 2, MPI_FLOAT, rightid, 0, MPI_COMM_WORLD);
                MPI_Recv(bufferL, size * 2, MPI_FLOAT, leftid, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(bufferR, size * 2, MPI_FLOAT, rightid, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                
                //read buffer
                for (int i = 1; i < size + 1; i++) {
                    for (int j = 0; j < 2; j++) {
                        arrayA[i * m + j] = bufferL[(i - 1) * 2 + j];
                    }
                }
                
                for (int i = 1; i < size + 1; i++) {
                    for (int j = 0; j < 2; j++) {
                        arrayA[i * m + j + blockSize + 2] = bufferR[(i - 1) * 2 + j];
                    }
                }
            }
        }

        
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
