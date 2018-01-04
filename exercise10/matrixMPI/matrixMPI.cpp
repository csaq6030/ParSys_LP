#include <iostream>
#include "mpi.h"
#define N 12


using namespace std;

int main(int argc, char **argv) {
    int numtasks,taskid,numworkers,source,dest,rows,offset,i,j,k;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Status status;
    
    numworkers = numtasks-1;
    
    
    rows = N / numtasks;
    
    //master process
    if (taskid == 0) {
        
        double* a = new double [N * N];
        double* b = new double [N * N];
        double* c = new double [N * N];
        
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                a[i * N + j] = 1.0;
                b[i * N + j] = 2.0;
                c[i * N + j] = -1.0; // for test purpose
            }
        }
        
        //send data
        offset = rows;
        
        for (dest = 1; dest <= numworkers; dest++) {
            MPI_Send(&a[offset * N], rows * N, MPI_DOUBLE,dest,1, MPI_COMM_WORLD);
            MPI_Send(b, N * N, MPI_DOUBLE, dest, 1, MPI_COMM_WORLD);
            offset += rows;
        }
        
        //own calculation
        
        for (k = 0; k < N; k++) {
            for (i = 0; i < rows; i++) {
                c[i * N + k] = 0.0;
                for (j = 0; j < N; j++)
                    c[i * N + k] += a[i * N + j] * b[j * N + k];
            }
        }
        
        //get data from other
        offset = rows;
        
        for (i = 1; i <= numworkers; i++) {
            source = i;
            MPI_Recv(&c[offset * N], rows * N, MPI_DOUBLE, source, 2, MPI_COMM_WORLD, &status);
            offset += rows;
        }
        
        printf("Result:\n");
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++)
                cout << c[i * N + j] << " ";
            cout << endl;
        }
        
        delete[] a;
        delete[] b;
        delete[] c;
        
    } else if (taskid > 0) { //other processes
        
        double* a = new double [rows * N];
        double* b = new double [N * N];
        double* c = new double [rows * N];
        
        source = 0;
        MPI_Recv(a, rows  * N, MPI_DOUBLE, source, 1, MPI_COMM_WORLD, &status);
        MPI_Recv(b, N * N, MPI_DOUBLE, source, 1, MPI_COMM_WORLD, &status);

        for (k = 0; k < N; k++) {
            for (i = 0; i < rows; i++) {
                c[i * N + k] = 0.0; // init to 0
                for (j = 0; j < N; j++)
                    c[i * N + k] += a[i * N + j] * b[j * N + k];
            }
        }
        
        MPI_Send(c, rows * N, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);
        
        delete[] a;
        delete[] b;
        delete[] c;
    }
    
    MPI_Finalize();
}
