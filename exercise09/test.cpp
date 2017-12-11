#include <iostream>
#include <cstdlib>
#include <mpi.h>
#include <cmath>
#include <algorithm>

using namespace std;

void printError() {
    cout << "Usage: [output] [size] [epsilon] [north] [south] [east] [west]" << endl;
}

void printArray2D(double *arr, int blockSize){
    for(int i = 0; i < blockSize; i++){
        for(int j = 0; j < blockSize; j++){
            cout << arr[i*blockSize + j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}


int main(int argc, char* argv[]) {
    int output = 0;
    int size = 0;
    double epsilon = 0.;
    double north = 0, south = 0, east = 0, west = 0;
    
    if (argc != 8) {
        printError();
        return EXIT_FAILURE;
    }
    output = atoi(argv[1]);
    size = atoi(argv[2]);
    epsilon = atof(argv[3]);
    north = atof(argv[4]);
    south = atof(argv[5]);
    east = atof(argv[6]);
    west = atof(argv[7]);
    
    int myid = 0, worldSize = 4;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
    // neighbor indexing for x*x=worldSize and initialize array
    //double *arrayA, *arrayB;
    MPI_Request req, rreq;
    if(myid == 0){
        double h = 5;
        MPI_Isend(&h,1,MPI_DOUBLE,1,MPI_ANY_TAG,MPI_COMM_WORLD,&req);  //west border
    }
    double d;
    if(myid !=0){  //top-left case
        double h =5;
        
        MPI_Irecv(&d,1,MPI_DOUBLE,0,MPI_ANY_TAG,MPI_COMM_WORLD,&rreq);
        
    }
    // calculate borders
    if(myid==0)
        MPI_Wait( &req, MPI_STATUSES_IGNORE);
    else
        MPI_Wait(&rreq, MPI_STATUSES_IGNORE);
    cout << myid << " " << d <<endl;

    MPI_Finalize();

    return EXIT_SUCCESS;
}
