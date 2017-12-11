#include <iostream>
#include <string>
#include <chrono>
#include <mpi.h>
#include <cmath>
#include <algorithm>

using namespace std;

void printError() {
    cout << "Usage: [output] [size] [epsilon] [north] [south] [east] [west]" << endl;
}

void printArray2D(double *arr, int blockSize, int id){
    for(int i = 0; i < blockSize; i++){
        cout << "id: " << id << "    " ;
        for(int j = 0; j < blockSize; j++){
            cout << arr[i*blockSize + j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}


int main(int argc, char* argv[]) {
    int iter = 0, output = 0, size, myid, worldSize;
    double epsilon = 0, progress = 0, reducedProgress = 0, north = 0, south = 0, east = 0, west = 0;
    
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
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
    
    // neighbor indexing for x*x=worldSize and initialize array
    int leftid, rightid, upid, downid;
    const int ij_range = static_cast<int>(sqrt(worldSize));
    const int ghostcells = 2;
    const int blockSize = (size/ij_range) + (ghostcells * 2);
    double *arrayA = NULL, *arrayB = NULL;
    
    if(myid == 0){  //top-left case
        leftid = MPI_PROC_NULL;
        upid = MPI_PROC_NULL;
        downid = myid + ij_range;
        rightid = myid + 1;
        arrayA = new double[blockSize * blockSize]();
        arrayB = new double[blockSize * blockSize]();
        for(int i = 0; i < blockSize; i++){
            arrayA[(blockSize) + i] = north;
            arrayA[i*(blockSize) + 1] = west;
            arrayB[(blockSize) + i] = north;
            arrayB[i*(blockSize) + 1] = west;
        }        
    }else if(myid < ij_range-1){ //top-middle cases
        upid = MPI_PROC_NULL;
        downid = myid + ij_range;
        leftid = myid - 1;
        rightid = myid + 1;
        
        arrayA = new double[blockSize * blockSize]();
        arrayB = new double[blockSize * blockSize]();
        for(int i = 0; i < blockSize; i++){
            arrayA[(blockSize) + i] = north;
            arrayB[(blockSize) + i] = north;
        }
    }else if(myid == ij_range-1){ // top-right case
        upid = MPI_PROC_NULL;
        rightid = MPI_PROC_NULL;
        downid = myid + ij_range;
        leftid = myid - 1;
        
        arrayA = new double[blockSize * blockSize]();
        arrayB = new double[blockSize * blockSize]();
        for(int i = 0; i < blockSize; i++){
            arrayA[(blockSize) + i] = north;
            arrayA[i*(blockSize)-2] = east;
            arrayB[(blockSize) + i] = north;
            arrayB[i*(blockSize)-2] = east;
        }
    }else if(myid == worldSize - ij_range){ // left-bottom case
        downid = MPI_PROC_NULL;
        leftid = MPI_PROC_NULL;
        upid = myid - ij_range;
        rightid = myid + 1;
        
        arrayA = new double[blockSize * blockSize]();
        arrayB = new double[blockSize * blockSize]();
        for(int i = 0; i < blockSize; i++){
            arrayA[i*(blockSize) + 1] = west;
            arrayA[blockSize*blockSize - i -1 -blockSize] = south;
            arrayB[i*(blockSize) + 1] = west;
            arrayB[blockSize*blockSize - i -1 -blockSize] = south;
        }
    }else if(myid == worldSize-1){    // right-bottom case
        rightid = MPI_PROC_NULL;
        downid = MPI_PROC_NULL;
        upid = myid - ij_range;
        leftid = myid - 1;
        
        arrayA = new double[blockSize * blockSize]();
        arrayB = new double[blockSize * blockSize]();
        for(int i = 0; i < blockSize; i++){
            arrayA[i*(blockSize)-2] = east;
            arrayA[blockSize*blockSize - i -1 -blockSize] = south;
            arrayB[i*(blockSize)-2] = east;
            arrayB[blockSize*blockSize - i -1 -blockSize] = south;
        }
    }else if(myid % ij_range == 0){ // left-middle cases
        leftid = MPI_PROC_NULL;
        upid = myid - ij_range;
        downid = myid + ij_range;
        rightid = myid + 1;
        
        arrayA = new double[blockSize * blockSize]();
        arrayB = new double[blockSize * blockSize]();
        for(int i = 0; i < blockSize; i++){
            arrayA[i*(blockSize) + 1] = west;
            arrayB[i*(blockSize) + 1] = west;
        }
    }else if(myid % ij_range == ij_range-1){  // right-middle cases
        rightid = MPI_PROC_NULL;
        upid = myid - ij_range;
        downid = myid + ij_range;
        leftid = myid - 1;
        
        arrayA = new double[blockSize * blockSize]();
        arrayB = new double[blockSize * blockSize]();
        for(int i = 0; i < blockSize; i++){
            arrayA[i*(blockSize)-2] = east;
            arrayB[i*(blockSize)-2] = east;
        }
    }else if(myid > ij_range*ij_range - ij_range){   // bottom-middle cases
        downid = MPI_PROC_NULL;
        upid = myid - ij_range;
        leftid = myid - 1;
        rightid = myid + 1;
        
        arrayA = new double[blockSize * blockSize]();
        arrayB = new double[blockSize * blockSize]();
        for(int i = 0; i < blockSize; i++){
            arrayA[blockSize*blockSize - i -1 -blockSize] = south;
            arrayB[blockSize*blockSize - i -1 -blockSize] = south;
        }
    }else{  //interiors
        upid = myid - ij_range;
        downid = myid + ij_range;
        leftid = myid - 1;
        rightid = myid + 1;
        
        arrayA = new double[blockSize * blockSize]();
        arrayB = new double[blockSize * blockSize]();
    }
    
    double *upBorder = new double[2*blockSize]();
    double *downBorder = new double[2*blockSize]();
    double *upRecvBuff = new double[2*blockSize]();
    double *downRecvBuff = new double[2*blockSize]();
    double *leftBorder = new double[2*(blockSize-4-4)]();
    double *rightBorder = new double[2*(blockSize-4-4)]();
    double *leftRecvBuff = new double[2*(blockSize-4-4)]();
    double *rightRecvBuff = new double[2*(blockSize-4-4)]();
    
    do{
        progress = 0;
        // calculate border values
        // first iteration, calculate only when neighbor exist -> ignore when borderline
        for(int i = 1; i < 5; i++){
            for(int j = 2; j < blockSize-2; j++){
                // top border cells
                if(i > 1 || (upid != MPI_PROC_NULL && i == 1)){
                    arrayB[i*blockSize + j] = (arrayA[i*blockSize + j] + arrayA[(i-1)*blockSize + j] + arrayA[(i+1)*blockSize + j] 
                                            + arrayA[i*blockSize + j - 1] + arrayA[i*blockSize + j + 1])*0.2;
                }
                // bottom border cells
                if(i > 1 || (downid != MPI_PROC_NULL && i == 1)){
                    arrayB[blockSize*blockSize - (i+1)*blockSize + j] = (arrayA[blockSize*blockSize - (i+1)*blockSize + j]
                                            + arrayA[blockSize*blockSize - (i+2)*blockSize + j] + arrayA[blockSize*blockSize - (i+1)*blockSize + j + 1] 
                                            + arrayA[blockSize*blockSize - (i)*blockSize + j] + arrayA[blockSize*blockSize - (i+1)*blockSize + j - 1])*0.2;
                }
            }
        }
        for(int j = 5; j < blockSize-4; j++){
            for(int i = 1; i < 5; i++){
                // left border cells
                if(i > 1 || (leftid != MPI_PROC_NULL && i == 1)){
                    arrayB[j*blockSize + i] = (arrayA[j*blockSize + i] + arrayA[(j-1)*blockSize + i] + arrayA[(j-1)*blockSize + i] 
                                            + arrayA[j*blockSize + i - 1] + arrayA[j*blockSize + i + 1])*0.2;
                }
                // right border cells
                if(i > 1 || (rightid != MPI_PROC_NULL && i == 1)){
                    arrayB[(j+1)*blockSize - 1 - i] = (arrayA[(j+1)*blockSize - 1 - i] + arrayA[(j)*blockSize - 1 - i] + arrayA[(j+2)*blockSize - 1 - i] 
                                            + arrayA[(j+1)*blockSize - 2 - i] + arrayA[(j+1)*blockSize - i])*0.2;
                }
            }
        }
        // 2nd iteration
        for(int i = 2; i < 4; i++){
            for(int j = 2; j < blockSize-2; j++){
                // top border cells
                arrayA[i*blockSize + j] = (arrayB[i*blockSize + j] + arrayB[(i-1)*blockSize + j] + arrayB[(i+1)*blockSize + j] 
                                            + arrayB[i*blockSize + j - 1] + arrayB[i*blockSize + j + 1])*0.2;
                // bottom border cells
                arrayA[blockSize*blockSize - (i+1)*blockSize + j] = (arrayB[blockSize*blockSize - (i+1)*blockSize + j]
                                            + arrayB[blockSize*blockSize - (i+2)*blockSize + j] + arrayB[blockSize*blockSize - (i+1)*blockSize + j + 1] 
                                            + arrayB[blockSize*blockSize - (i)*blockSize + j] + arrayB[blockSize*blockSize - (i+1)*blockSize + j - 1])*0.2;
                progress += std::abs(arrayB[i*blockSize + j] - arrayA[i*blockSize + j]);
                progress += std::abs(arrayB[blockSize*blockSize - (i+1)*blockSize + j] - arrayA[blockSize*blockSize - (i+1)*blockSize + j]);
            }
        }
        for(int j = 4; j < blockSize-4; j++){
            for(int i = 2; i < 4; i++){
                // left border cells
                arrayA[j*blockSize + i] = (arrayB[j*blockSize + i] + arrayB[(j-1)*blockSize + i] + arrayB[(j-1)*blockSize + i] 
                                            + arrayB[j*blockSize + i - 1] + arrayB[j*blockSize + i + 1])*0.2;
                // right border cells
                arrayA[(j+1)*blockSize - 1 - i] = (arrayB[(j+1)*blockSize - 1 - i] + arrayB[(j)*blockSize - 1 - i] + arrayB[(j+2)*blockSize - 1 - i] 
                                            + arrayB[(j+1)*blockSize - 2 - i] + arrayB[(j+1)*blockSize - i])*0.2;
                progress += std::abs(arrayB[j*blockSize + i] - arrayA[j*blockSize + i]);
                progress += std::abs(arrayB[(j+1)*blockSize - 1 - i] - arrayA[(j+1)*blockSize - 1 - i]);
            }
        }
        // create send buffers
        memcpy(upBorder, &arrayA[2*blockSize], 2*blockSize*sizeof(double));
        memcpy(downBorder, &arrayA[blockSize*blockSize - 4*blockSize], 2*blockSize*sizeof(double));
        
        for(int i = 0; i < blockSize-4-4; i++){
            memcpy(&leftBorder[i*2], &arrayA[(i+4)*blockSize+2], 2*sizeof(double));
            memcpy(&rightBorder[i*2], &arrayA[(i+5)*blockSize-4], 2*sizeof(double));
        }
        // send out/receive ghost-cells
        MPI_Request req[8];
        MPI_Isend(upBorder,2*blockSize,MPI_DOUBLE,upid,0,MPI_COMM_WORLD,&req[0]); // north border
        MPI_Isend(downBorder,2*blockSize,MPI_DOUBLE,downid,0,MPI_COMM_WORLD,&req[1]);  //south border
        MPI_Isend(rightBorder,2*(blockSize-4-4),MPI_DOUBLE,rightid,0,MPI_COMM_WORLD,&req[2]);  //east border
        MPI_Isend(leftBorder,2*(blockSize-4-4),MPI_DOUBLE,leftid,0,MPI_COMM_WORLD,&req[3]);  //west border

        MPI_Irecv(leftRecvBuff,2*(blockSize-4-4),MPI_DOUBLE,leftid,0,MPI_COMM_WORLD,&req[6]);
        MPI_Irecv(rightRecvBuff,2*(blockSize-4-4),MPI_DOUBLE,rightid,0,MPI_COMM_WORLD,&req[7]);
        MPI_Irecv(upRecvBuff,2*blockSize,MPI_DOUBLE,upid,0,MPI_COMM_WORLD,&req[4]);
        MPI_Irecv(downRecvBuff,2*blockSize,MPI_DOUBLE,downid,0,MPI_COMM_WORLD,&req[5]);

        // calculate interiors
        for(int i = 5; i < blockSize-5; i++){
            for(int j = 5; j < blockSize-5; j++){
                arrayB[i*blockSize + j] = (arrayA[i*blockSize + j] + arrayA[(i-1)*blockSize + j] + arrayA[(i+1)*blockSize + j] 
                                            + arrayA[i*blockSize + j - 1] + arrayA[i*blockSize + j + 1])*0.2;
            }
        }
        for(int i = 4; i < blockSize-4; i++){
            for(int j = 4; j < blockSize-4; j++){
                arrayA[i*blockSize + j] = (arrayB[i*blockSize + j] + arrayB[(i-1)*blockSize + j] + arrayB[(i+1)*blockSize + j] 
                                            + arrayB[i*blockSize + j - 1] + arrayB[i*blockSize + j + 1])*0.2;
                progress += std::abs(arrayA[i*blockSize + j] - arrayB[i*blockSize + j]);
            }
        }

        MPI_Allreduce(&progress, &reducedProgress, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        // load left and right borders -> TODO use MPI_Test()
        MPI_Wait(&req[6], MPI_STATUSES_IGNORE);
        MPI_Wait(&req[7], MPI_STATUSES_IGNORE);
        // TODO can be avoided with MPI type
        for(int i = 0; i < blockSize-4-4; i++){
            if(leftid!=MPI_PROC_NULL)
                memcpy(&arrayA[(i+4)*blockSize], &leftRecvBuff[i*2], 2*sizeof(double));
            if(rightid!=MPI_PROC_NULL)
                memcpy(&arrayA[(i+5)*blockSize-2], &rightRecvBuff[i*2], 2*sizeof(double));
        }
        
//        for(int i = 0; i <= blockSize-4-4; i++){
//            if(leftid!=MPI_PROC_NULL){
//                arrayA[(i+4)*blockSize]=leftRecvBuff[i*2];
//                arrayA[(i+4)*blockSize+1]=leftRecvBuff[i*2+1];
//            }
//            if(rightid!=MPI_PROC_NULL){
//                arrayA[(i+5)*blockSize-2]=rightRecvBuff[i*2];
//                arrayA[(i+5)*blockSize-2+1]=rightRecvBuff[i*2+1];
//            }
//        }
        
//        for (int i = 0; i < blockSize - 8; i++) {
//            for (int j = 0; i < 2; j++) {
//                arrayA[(i+4)*blockSize] = rightRecvBuff[i * 2 + j];
//            }
//        }
        MPI_Wait(&req[4], MPI_STATUSES_IGNORE);
        if(upid!=MPI_PROC_NULL)
            memcpy(&arrayA[0], upRecvBuff, 2*blockSize*sizeof(double));
        MPI_Wait(&req[5], MPI_STATUSES_IGNORE);        
        if(downid!=MPI_PROC_NULL)
            memcpy(&arrayA[blockSize*blockSize-2*blockSize], downRecvBuff, 2*blockSize*sizeof(double));
        
        if(myid == 0){
            printArray2D(rightRecvBuff,2*(blockSize-8), 0);
        }
        if(myid==0)
            printArray2D(arrayA, blockSize, myid);

        // synchronize 
        MPI_Waitall(4, req, MPI_STATUSES_IGNORE);
        
        if(myid==0)
            iter +=2;
    }while(reducedProgress >= epsilon);
    
    MPI_Finalize();
    
    if(myid==0){
        cout << "Iterations: " << iter << endl;
    }
    cout << myid <<endl;
    
    
    delete[] arrayA;
    delete[] arrayB;
    delete[] upBorder;
    delete[] downBorder;
    delete[] leftBorder;
    delete[] rightBorder;
    delete[] upRecvBuff;
    delete[] downRecvBuff;
    delete[] leftRecvBuff;
    delete[] rightRecvBuff;

    return EXIT_SUCCESS;
}
