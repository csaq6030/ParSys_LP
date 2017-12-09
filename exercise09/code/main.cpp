#include <iostream>
#include <cmath>
#include <mpi.h>

using namespace std;


inline void border1Array(double **arrayA, double **arrayB,  int size, const double up, double down, const double  left, const double right) {
    for (int i = 0; i < size; i++) {
        arrayA[0][i] = up;
        arrayA[i][0] = left;
        arrayA[size - 1][i] = down;
        arrayA[i][size - 1] = right;
        
        arrayB[0][i] = up;
        arrayB[i][0] = left;
        arrayB[size - 1][i] = down;
        arrayB[i][size - 1] = right;
    }
}

inline void print2Darray(double **array, const int start_i, const int start_j , const int n, const int m) {
    for (int i = start_i; i < n; i++) {
        for (int j = start_j; j < m; j++) {
            cout << array[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}

inline double stencil1Array(double **arrayA, double **arrayB, const int size) {
    double progress = 0;
    
    for (int i = 1;  i < size - 1; i++) {
        for (int j = 1; j < size - 1; j++) {
            arrayB[i][j] = (arrayA[i - 1][j] + arrayA[i + 1][j]
                            + arrayA[i][j - 1] + arrayA[i][j] + arrayA[i][j + 1]) / 5;
        }
    }
    
    for (int i = 1;  i < size - 1; i++) {
        for (int j = 1; j < size - 1; j++) {
            arrayA[i][j] = (arrayB[i - 1][j] + arrayB[i + 1][j]
                            + arrayB[i][j - 1] + arrayB[i][j] + arrayB[i][j + 1]) / 5;
            progress += std::abs(arrayA[i][j] - arrayB[i][j]);
        }
    }
    
    return progress;
    
}

inline void print2DarrayBlock (double *array, const int start_i, const int start_j , const int n, const int m) {
    for (int i = start_i; i < n; i++) {
        for (int j = start_j; j < m; j++) {
            cout << array[i * m + j] << " ";
        }
        cout << endl;
    }
    cout << endl;
    
}

inline double stencil2D(double *arrayB, double *arrayA, const int m, const int start_i_1, const int start_j_1,
                 const int end_i_1, const int end_j_1, const int start_i_2, const int start_j_2, const int end_i_2,
                 const int end_j_2) {
    double progress = 0;
    for (int i = start_i_1; i < end_i_1; i++){
                    for (int j = start_j_1; j < end_j_1; j++) {
                        arrayB[i * m + j] = (arrayA[i * m + j - 1] + arrayA[i * m + j] + arrayA[i * m + j + 1] +
                                             arrayA[(i - 1) * m + j] + arrayA[(i + 1) * m + j])/5;
                    }
                }

    for (int i = start_i_2; i < end_i_2; i++){
                    for (int j = start_j_2; j < end_j_2; j++) {
                        arrayA[i * m + j] = (arrayB[i * m + j - 1] + arrayB[i * m + j] + arrayB[i * m + j + 1] +
                                             arrayB[(i - 1) * m + j] + arrayB[(i + 1) * m + j])/5;
                        progress += std::abs(arrayB[i * m + j] - arrayA[i * m + j]);
                    }
                }
    return progress;
}

int main(int argc, char* argv[]) {
    const bool output = true;
    
    const int size = 10; //512 or 768
    
    const double up = 1.;
    const double down = 0.;
    const double left = -0.5;
    const double right = 0.5;
    double reducedProgress = 0.;
    
    const double epsilon = 1.; // 10 or 100 for 512 or 768 according
    
    //init mpi
    int myid, worldSize;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
    
    if (worldSize == 1) {
        
        double **arrayA=new double*[size + 2];
        double **arrayB=new double*[size + 2];
        for( int i = 0; i< size + 2; ++i ) { //i
            arrayA[i] = new double[size + 2];
            arrayB[i] = new double[size + 2];
            for( int j=0; j< size + 2; ++j ) {//j
                arrayA[i][j] = 0;
                arrayB[i][j] = 0;
            }
        }
        
        border1Array(arrayA, arrayB, size + 2, up, down, left, right);
        //print2Darray(arrayA, 0, 0,size + 2, size + 2);
        
        while (stencil1Array(arrayA, arrayB, size + 2) >= epsilon) {
            }
        
        if (output) {
            print2Darray(arrayA, 1, 1,size + 1, size + 1);
        }
        
    
    
        for(int i = 0; i < size + 2; i++) {
            delete[] arrayA[i];
            delete[] arrayB[i];
        }
        delete[] arrayA;
        delete[] arrayB;
        
        
    } else if (worldSize == 2) {
        const int blockSize = size / 2;
        const int n = size + 2; //m equals i in A[i][j]
        const int m = (blockSize + 3); //m equals j in A[i][j]
        
        if (myid == 0) { //left
            
            double *arrayA = new double[n * m]();
            double *arrayB = new double[n * m]();
            //boundary left
            for (int i = 0; i < n; i++) {
                arrayA[i * m] = left;
                arrayB[i * m] = left;
            }
            
            //boundary up and down
            for (int i = 0; i < m; i++) {
                arrayA[i] = up;
                arrayB[i] = up;
                
                arrayA[(n-1) * m + i] = down;
                arrayB[(n-1) * m + i] = down;
            }
            
            //cout << "done init first" << endl;
            double *buffer = new double[size * 2]();
            
            while(true){

                double progress = stencil2D(arrayB, arrayA, m, 1, 1, n - 1, m - 1, 1,
                                            1, n - 1, m - 2);
                //cout << "iter first" << endl;
                MPI_Allreduce(&progress, &reducedProgress, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                if(reducedProgress < epsilon){
                    //send result - none master //collect data
                    //own
                    double *result = new double[size * size]();
                    for(int i = 1; i < size + 1; i++){
                        for (int j = 1; j < blockSize + 1; j++) {
                            result[(i - 1) * size + (j - 1)] = arrayA[i * m + j];
                        }
                    }
                    //other
                    
                    double *tmp2 = new double[size * blockSize]();
                    //last result
                    MPI_Recv(tmp2,size * blockSize, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    for(int i = 0; i < size ; i++){
                        for (int j = 0; j < blockSize; j++) {
                            result[i * size + j + blockSize] = tmp2[i * blockSize + j];
                            //result[i * size + j + (worldSize - 1) * blockSize] = 1111; //little test
                        }
                    }
                    delete[] tmp2;
                    
                    //Print result
                    if (output) {
                        print2DarrayBlock(result, 0, 0, size, size);
                    }
                    
                    delete[] buffer;
                    delete[] result;
                    delete[] arrayA;
                    delete[] arrayB;
                    break;
                } else {
                
                    //fill buffer
                    for (int i = 1; i < size + 1; i++) {
                        for (int j = blockSize - 1; j < blockSize + 1; j++) {
                            buffer[(i - 1) * 2 + j - blockSize + 1] = arrayA[i * m + j];
                        }
                    }
                    
                    MPI_Send(buffer, size * 2, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
                    MPI_Recv(buffer, size * 2, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    
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
            }
                
                
        } else { // right
            double *arrayA = new double[n * m]();
            double *arrayB = new double[n * m]();
            
            
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
            
            double *buffer = new double[size * 2]();
            
            //cout << "done init last" << endl;
            
            while(true){

                double progress = stencil2D(arrayB, arrayA, m, 1, 1, n - 1, m - 1, 1,
                                            2, n - 1, m - 1);

                //cout << "iter last" << endl;
                MPI_Allreduce(&progress, &reducedProgress, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                if(reducedProgress < epsilon){
                    //cout << "result last sending" << endl;
                    //send result
                    double *tmp = new double[size * blockSize]();
                    
                    for (int i = 1; i < size + 1; i++) {
                        for (int j = 2; j < m - 1; j++) {
                            tmp[(i - 1) * blockSize + (j - 2)] = arrayA[i * m + j];
                        }
                    }
                    
                    
                    
                    MPI_Send(tmp, blockSize * size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
                    
                    delete[] tmp;
                    delete[] buffer;
                    delete[] arrayA;
                    delete[] arrayB;
                    break;
                } else {
                    //cout << "fill buffer last" << endl;
                    //fill buffer
                    for (int i = 1; i < size + 1; i++) {
                        for (int j = 2; j < 4; j++) {
                            buffer[(i - 1) * 2 + j - 2] = arrayA[i * m + j];
                        }
                    }
                    //cout << "done fill buffer last, now sending" << endl;
                    MPI_Send(buffer, size * 2, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
                    //cout << "send last, now recv" << endl;
                    MPI_Recv(buffer, size * 2, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    //cout << "done recv last, now reading buffer" << endl;
                    
                    //read buffer
                    for (int i = 1; i < size + 1; i++) {
                        for (int j = 0; j < 2; j++) {
                            //cout << "i: " << i << " " << "j: " << j << endl;
                            arrayA[i * m + j] = buffer[(i - 1) * 2 + j];
                        }
                    }
                    
                    //cout << "done reading buffer last" << endl;
                }
            }
            
        }
    } else { // othe worldsizes 2 cases missing square ones (4,16,64) and m != n ones (8 (4*2 || 2* 4), 32 (8*4 || 4 * 8))


    }
    
    
    
    
    MPI_Finalize();
    return EXIT_SUCCESS;
    
}
