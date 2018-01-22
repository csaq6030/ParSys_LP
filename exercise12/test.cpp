#include <iostream>
#include "mpi.h"
#include <omp.h>


using namespace std;

int main(int argc, char **argv) {
    
    int numtasks,taskid;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Status status;
    
    
#pragma openmp parallel
    {
        int threadNum = omp_get_thread_num();
        cout << "mpirank: " << numtasks << "num omp thread: " << thread << endl;
    }
    
    
    MPI_Finalize();
    
    
}
