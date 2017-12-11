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

inline void print2Darray(double **array, const int start_i, const int start_j , const int n, const int m, const int id) {
    for (int i = start_i; i < n; i++) {
        cout << "id: " << id << "    " ;
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

inline void print2Darray(double *array, const int start_i, const int start_j , const int n, const int m, const int id) {
    for (int i = start_i; i < n; i++) {
        cout << "id: " << id << "    " ;
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
    
    const int size = 16; //512 or 768
    
    const double up = 1.;
    const double down = 0.;
    const double left = -0.5;
    const double right = 0.5;
    double reducedProgress = 0.;
    int iter = 0;
    
    const double epsilon = 100.; // 10 or 100 for 512 or 768 according
    
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

            iter += 2;

            if (iter == 14)
                break;
        }
        
        if (output) {
            print2Darray(arrayA, 1, 1,size + 1, size + 1, 0);

            cout << "iter: " << iter << endl;
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
                        print2Darray(result, 0, 0, size, size, 0);
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
    } else if (int (sqrt(worldSize)) * int (sqrt(worldSize)) == worldSize){ //square case
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
                arrayA[(blockSize) + i] = up;
                arrayA[i*(blockSize) + 1] = left;
                arrayB[(blockSize) + i] = up;
                arrayB[i*(blockSize) + 1] = left;
            }
        }else if(myid < ij_range-1){ //top-middle cases
            upid = MPI_PROC_NULL;
            downid = myid + ij_range;
            leftid = myid - 1;
            rightid = myid + 1;

            arrayA = new double[blockSize * blockSize]();
            arrayB = new double[blockSize * blockSize]();
            for(int i = 0; i < blockSize; i++){
                arrayA[(blockSize) + i] = up;
                arrayB[(blockSize) + i] = up;
            }
        }else if(myid == ij_range-1){ // top-right case
            upid = MPI_PROC_NULL;
            rightid = MPI_PROC_NULL;
            downid = myid + ij_range;
            leftid = myid - 1;

            arrayA = new double[blockSize * blockSize]();
            arrayB = new double[blockSize * blockSize]();
            for(int i = 0; i < blockSize; i++){
                arrayA[(blockSize) + i] = up;
                arrayA[i*(blockSize)-2] = right;
                arrayB[(blockSize) + i] = up;
                arrayB[i*(blockSize)-2] = right;
            }
        }else if(myid == worldSize - ij_range){ // left-bottom case
            downid = MPI_PROC_NULL;
            leftid = MPI_PROC_NULL;
            upid = myid - ij_range;
            rightid = myid + 1;

            arrayA = new double[blockSize * blockSize]();
            arrayB = new double[blockSize * blockSize]();
            for(int i = 0; i < blockSize; i++){
                arrayA[i*(blockSize) + 1] = left;
                arrayA[blockSize*blockSize - i -1 -blockSize] = down;
                arrayB[i*(blockSize) + 1] = left;
                arrayB[blockSize*blockSize - i -1 -blockSize] = down;
            }
        }else if(myid == worldSize-1){    // right-bottom case
            rightid = MPI_PROC_NULL;
            downid = MPI_PROC_NULL;
            upid = myid - ij_range;
            leftid = myid - 1;

            arrayA = new double[blockSize * blockSize]();
            arrayB = new double[blockSize * blockSize]();
            for(int i = 0; i < blockSize; i++){
                arrayA[i*(blockSize)-2] = right;
                arrayA[blockSize*blockSize - i -1 -blockSize] = down;
                arrayB[i*(blockSize)-2] = right;
                arrayB[blockSize*blockSize - i -1 -blockSize] = down;
            }
        }else if(myid % ij_range == 0){ // left-middle cases
            leftid = MPI_PROC_NULL;
            upid = myid - ij_range;
            downid = myid + ij_range;
            rightid = myid + 1;

            arrayA = new double[blockSize * blockSize]();
            arrayB = new double[blockSize * blockSize]();
            for(int i = 0; i < blockSize; i++){
                arrayA[i*(blockSize) + 1] = left;
                arrayB[i*(blockSize) + 1] = left;
            }
        }else if(myid % ij_range == ij_range-1){  // right-middle cases
            rightid = MPI_PROC_NULL;
            upid = myid - ij_range;
            downid = myid + ij_range;
            leftid = myid - 1;

            arrayA = new double[blockSize * blockSize]();
            arrayB = new double[blockSize * blockSize]();
            for(int i = 0; i < blockSize; i++){
                arrayA[i*(blockSize)-2] = right;
                arrayB[i*(blockSize)-2] = right;
            }
        }else if(myid > ij_range*ij_range - ij_range){   // bottom-middle cases
            downid = MPI_PROC_NULL;
            upid = myid - ij_range;
            leftid = myid - 1;
            rightid = myid + 1;

            arrayA = new double[blockSize * blockSize]();
            arrayB = new double[blockSize * blockSize]();
            for(int i = 0; i < blockSize; i++){
                arrayA[blockSize*blockSize - i -1 -blockSize] = down;
                arrayB[blockSize*blockSize - i -1 -blockSize] = down;
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
        double *leftBorder = new double[2*(blockSize-4)]();
        double *rightBorder = new double[2*(blockSize-4)]();
        double *leftRecvBuff = new double[2*(blockSize-4)]();
        double *rightRecvBuff = new double[2*(blockSize-4)]();

        double progress;
        do{
            progress = 0;
            // calculate border values
            // first iteration, calculate only when neighbor exist -> ignore when borderline
            for(int i = 1; i < 5; i++){
                for(int j = 1; j < blockSize-1; j++){
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
            for(int j = 5; j < blockSize-5; j++){
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

            if (output && myid == 0) {
                print2Darray(arrayB, 0, 0, blockSize, blockSize, myid);
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

            for(int i = 0; i < blockSize-4; i++){
                memcpy(&leftBorder[i*2], &arrayA[(i+2)*blockSize+2], 2*sizeof(double));
                memcpy(&rightBorder[i*2], &arrayA[(i+3)*blockSize-4], 2*sizeof(double));
            }
            // send out/receive ghost-cells
            MPI_Request req[8];
            MPI_Isend(upBorder,2*blockSize,MPI_DOUBLE,upid,0,MPI_COMM_WORLD,&req[0]); // up border
            MPI_Isend(downBorder,2*blockSize,MPI_DOUBLE,downid,0,MPI_COMM_WORLD,&req[1]);  //down border
            MPI_Isend(rightBorder,2*(blockSize-4),MPI_DOUBLE,rightid,0,MPI_COMM_WORLD,&req[2]);  //right border
            MPI_Isend(leftBorder,2*(blockSize-4),MPI_DOUBLE,leftid,0,MPI_COMM_WORLD,&req[3]);  //left border

            MPI_Irecv(leftRecvBuff,2*(blockSize-4),MPI_DOUBLE,leftid,0,MPI_COMM_WORLD,&req[6]);
            MPI_Irecv(rightRecvBuff,2*(blockSize-4),MPI_DOUBLE,rightid,0,MPI_COMM_WORLD,&req[7]);
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
            for(int i = 0; i < blockSize-4; i++){
                if(leftid!=MPI_PROC_NULL)
                    memcpy(&arrayA[(i+2)*blockSize], &leftRecvBuff[i*2], 2*sizeof(double));
                if(rightid!=MPI_PROC_NULL)
                    memcpy(&arrayA[(i+3)*blockSize-2], &rightRecvBuff[i*2], 2*sizeof(double));
            }

            MPI_Wait(&req[4], MPI_STATUSES_IGNORE);
            if(upid!=MPI_PROC_NULL)
                memcpy(&arrayA[0], upRecvBuff, 2*blockSize*sizeof(double));
            MPI_Wait(&req[5], MPI_STATUSES_IGNORE);
            if(downid!=MPI_PROC_NULL)
                memcpy(&arrayA[blockSize*blockSize-2*blockSize], downRecvBuff, 2*blockSize*sizeof(double));

            if(myid==0)
                iter +=2;

            // synchronize 
            MPI_Waitall(4, req, MPI_STATUSES_IGNORE);
        }while(reducedProgress >= epsilon);



        if(myid==0 && output){
            cout << "Iterations: " << iter << endl;
        }
        //cout << myid <<endl;
        if (output && myid == 0) {
            print2Darray(arrayA, 0, 0, blockSize, blockSize, myid);
        }



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


    } else { //not square ones (8 (2*4), 32 (2*16))
        const int columnSize = worldSize / 2;
        int leftid, rightid, upid, downid;

        if (myid == 0) { //  upper left corner
            leftid = MPI_PROC_NULL;
            rightid = 1;
            upid = MPI_PROC_NULL;
            downid = columnSize;
        } else if ( myid < columnSize - 1){ // upper row minus right corner
            leftid = myid - 1;
            rightid = myid + 1;
            upid = MPI_PROC_NULL;
            downid = myid + columnSize;
        } else if (myid == columnSize - 1 ) { // upper right corner
            leftid = myid - 1;
            rightid = MPI_PROC_NULL;
            upid = MPI_PROC_NULL;
            downid = myid + columnSize;
        } else if (myid == columnSize) { // lower left corner
            leftid = MPI_PROC_NULL;
            rightid = myid + 1;
            upid = 0;
            downid = MPI_PROC_NULL;
        } else if (myid < worldSize - 1) { // lower row minus right corner
            leftid = myid - 1;
            rightid = myid + 1;
            upid = myid - columnSize;
            downid = MPI_PROC_NULL;
        } else { // right corner
            leftid = myid - 1;
            rightid = MPI_PROC_NULL;
            upid = myid - columnSize;
            downid = MPI_PROC_NULL;
        }



    }
    
    
    
    
    MPI_Finalize();
    return EXIT_SUCCESS;
    
}
