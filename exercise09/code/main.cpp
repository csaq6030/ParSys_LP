#include <iostream>
#include <cmath>
#include <chrono>
#include <mpi.h>
#include "helper.h"

using namespace std;

int main(int argc, char *argv[]) {
    const bool output = true;

    const int size = 512; //512 or 768

    const double up = 1.;
    const double down = 0;//0.;
    const double left = -0.5;//-0.5;
    const double right = 0.5;//0.5;
    double reducedProgress = 0.;
    int iter = 0;

    const double epsilon = 1.; // 10 or 100 for 512 or 768 according

    cout.precision(4);
    cout << fixed;

    int myid, worldSize;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);

    if (worldSize == 1) {

        double **arrayA = new double *[size + 2];
        double **arrayB = new double *[size + 2];
        for (int i = 0; i < size + 2; ++i) { //i
            arrayA[i] = new double[size + 2];
            arrayB[i] = new double[size + 2];
            for (int j = 0; j < size + 2; ++j) {//j
                arrayA[i][j] = 0;
                arrayB[i][j] = 0;
            }
        }

        border1Array(arrayA, arrayB, size + 2, up, down, left, right);
        //print2Darray(arrayA, 0, 0,size + 2, size + 2);

        // start time measurement for the stencil operation
        const std::chrono::time_point<std::chrono::high_resolution_clock>
        start(std::chrono::high_resolution_clock::now());
        do {
            iter += 2;
        } while (stencil1Array(arrayA, arrayB, size + 2) >= epsilon);

        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::high_resolution_clock::now() - start);
        cout << elapsed.count() << endl;

        if (output) {
            //print2Darray(arrayA, 0, 0,size + 2, size + 2, 0);
            cout << "iter: " << iter << endl;
        }

        for (int i = 0; i < size + 2; i++) {
            delete[] arrayA[i];
            delete[] arrayB[i];
        }
        delete[] arrayA;
        delete[] arrayB;

    } else if (worldSize % 2 == 1) {
        MPI_Finalize();
        if (myid == 0)
            cout << "odd # of processes not supported!!!" << endl;
        return EXIT_FAILURE;
    } else { // works for all cases where x*x=worldSize is fulfilled and 2,8,32 special cases
        int leftid, rightid, upid, downid;
        int i_range = 2;    //default case is for np = 8
        int j_range = 4;
        if (worldSize == 2) {
            i_range = 1;
            j_range = 2;
        } else if (worldSize == 32) {
            i_range = 4;
            j_range = 8;
        } else if (worldSize == int(sqrt(worldSize)) * int(sqrt(worldSize))) {
            i_range = int(sqrt(worldSize));
            j_range = i_range;
        }

        const int ghostcells = 2;
        const int iblockSize = (size / i_range) + (ghostcells * 2);
        const int jblockSize = (size / j_range) + (ghostcells * 2);
        double *arrayA = NULL, *arrayB = NULL;
        // MPI datatype to avoid buffers
        MPI_Datatype sideBorderType;
        MPI_Type_vector(jblockSize - 4, 2, iblockSize, MPI_DOUBLE, &sideBorderType);
        MPI_Type_commit(&sideBorderType);

        //neighbor indexing and initialize array
        if (myid == 0 && worldSize == 2) {   //special case with 3 borders
            leftid = MPI_PROC_NULL;
            upid = MPI_PROC_NULL;
            downid = myid + i_range;
            rightid = MPI_PROC_NULL;
            arrayA = new double[iblockSize * jblockSize]();
            arrayB = new double[iblockSize * jblockSize]();
            for (int i = 0; i < iblockSize; i++) {
                arrayA[(iblockSize) + i] = up;
                arrayB[(iblockSize) + i] = up;
                if (i < jblockSize) {
                    arrayA[i * (iblockSize) + 1] = left;
                    arrayB[i * (iblockSize) + 1] = left;
                    arrayA[(i + 1) * (iblockSize) - 2] = right;
                    arrayB[(i + 1) * (iblockSize) - 2] = right;
                }
            }
        } else if (myid == worldSize - 1 && worldSize == 2) { //special case with 3 borders
            rightid = MPI_PROC_NULL;
            downid = MPI_PROC_NULL;
            upid = myid - i_range;
            leftid = MPI_PROC_NULL;
            arrayA = new double[iblockSize * jblockSize]();
            arrayB = new double[iblockSize * jblockSize]();
            for (int i = 0; i < iblockSize; i++) {
                arrayA[iblockSize * jblockSize - i - 1 - iblockSize] = down;
                arrayB[iblockSize * jblockSize - i - 1 - iblockSize] = down;
                if (i < jblockSize) {
                    arrayB[(i + 1) * (iblockSize) - 2] = right;
                    arrayA[(i + 1) * (iblockSize) - 2] = right;
                    arrayA[i * (iblockSize) + 1] = left;
                    arrayB[i * (iblockSize) + 1] = left;
                }
            }
        } else if (myid == 0) {  //top-left case
            leftid = MPI_PROC_NULL;
            upid = MPI_PROC_NULL;
            downid = myid + i_range;
            rightid = myid + 1;
            arrayA = new double[iblockSize * jblockSize]();
            arrayB = new double[iblockSize * jblockSize]();
            for (int i = 0; i < iblockSize; i++) {
                arrayA[(iblockSize) + i] = up;
                arrayB[(iblockSize) + i] = up;
                if (i < jblockSize) {
                    arrayA[i * (iblockSize) + 1] = left;
                    arrayB[i * (iblockSize) + 1] = left;
                }
            }
        } else if (myid < i_range - 1) { //top-middle cases
            upid = MPI_PROC_NULL;
            downid = myid + i_range;
            leftid = myid - 1;
            rightid = myid + 1;
            arrayA = new double[iblockSize * jblockSize]();
            arrayB = new double[iblockSize * jblockSize]();
            for (int i = 0; i < iblockSize; i++) {
                arrayA[(iblockSize) + i] = up;
                arrayB[(iblockSize) + i] = up;
            }
        } else if (myid == i_range - 1) { // top-right case
            upid = MPI_PROC_NULL;
            rightid = MPI_PROC_NULL;
            downid = myid + i_range;
            leftid = myid - 1;

            arrayA = new double[iblockSize * jblockSize]();
            arrayB = new double[iblockSize * jblockSize]();
            for (int i = 0; i < iblockSize; i++) {
                arrayA[(iblockSize) + i] = up;
                arrayB[(iblockSize) + i] = up;
                if (i < jblockSize) {
                    arrayA[(i + 1) * (iblockSize) - 2] = right;
                    arrayB[(i + 1) * (iblockSize) - 2] = right;
                }
            }
        } else if (myid == worldSize - i_range) { // left-bottom case
            downid = MPI_PROC_NULL;
            leftid = MPI_PROC_NULL;
            upid = myid - i_range;
            rightid = myid + 1;

            arrayA = new double[iblockSize * jblockSize]();
            arrayB = new double[iblockSize * jblockSize]();
            for (int i = 0; i < iblockSize; i++) {
                arrayA[iblockSize * jblockSize - i - 1 - iblockSize] = down;
                arrayB[iblockSize * jblockSize - i - 1 - iblockSize] = down;
                if (i < jblockSize) {
                    arrayA[i * (iblockSize) + 1] = left;
                    arrayB[i * (iblockSize) + 1] = left;
                }
            }
        } else if (myid == worldSize - 1) {    // right-bottom case
            rightid = MPI_PROC_NULL;
            downid = MPI_PROC_NULL;
            upid = myid - i_range;
            leftid = myid - 1;
            arrayA = new double[iblockSize * jblockSize]();
            arrayB = new double[iblockSize * jblockSize]();
            for (int i = 0; i < iblockSize; i++) {
                arrayA[iblockSize * jblockSize - i - 1 - iblockSize] = down;
                arrayB[iblockSize * jblockSize - i - 1 - iblockSize] = down;
                if (i < jblockSize) {
                    arrayB[(i + 1) * (iblockSize) - 2] = right;
                    arrayA[(i + 1) * (iblockSize) - 2] = right;
                }
            }
        } else if (myid % i_range == 0) { // left-middle cases
            leftid = MPI_PROC_NULL;
            upid = myid - i_range;
            downid = myid + i_range;
            rightid = myid + 1;

            arrayA = new double[iblockSize * jblockSize]();
            arrayB = new double[iblockSize * jblockSize]();
            for (int i = 0; i < jblockSize; i++) {
                arrayA[i * (iblockSize) + 1] = left;
                arrayB[i * (iblockSize) + 1] = left;
            }
        } else if (myid % i_range == i_range - 1) {  // right-middle cases
            rightid = MPI_PROC_NULL;
            upid = myid - i_range;
            downid = myid + i_range;
            leftid = myid - 1;

            arrayA = new double[iblockSize * jblockSize]();
            arrayB = new double[iblockSize * jblockSize]();
            for (int i = 0; i < jblockSize; i++) {
                arrayA[(i + 1) * (iblockSize) - 2] = right;
                arrayB[(i + 1) * (iblockSize) - 2] = right;
            }
        } else if (myid > i_range * j_range - i_range) {   // bottom-middle cases
            downid = MPI_PROC_NULL;
            upid = myid - i_range;
            leftid = myid - 1;
            rightid = myid + 1;

            arrayA = new double[iblockSize * jblockSize]();
            arrayB = new double[iblockSize * jblockSize]();
            for (int i = 0; i < iblockSize; i++) {
                arrayA[iblockSize * jblockSize - i - 1 - iblockSize] = down;
                arrayB[iblockSize * jblockSize - i - 1 - iblockSize] = down;
            }
        } else {  //interiors
            upid = myid - i_range;
            downid = myid + i_range;
            leftid = myid - 1;
            rightid = myid + 1;

            arrayA = new double[iblockSize * jblockSize]();
            arrayB = new double[iblockSize * jblockSize]();
        }

		MPI_Barrier(MPI_COMM_WORLD);  
        const std::chrono::time_point<std::chrono::high_resolution_clock>
        start(std::chrono::high_resolution_clock::now());
        double progress;
        do {
            progress = 0;
            // calculate border values
            // first iteration, calculate only when neighbor exist -> ignore when borderline
            for (int i = 1; i < 5; i++) {
                int j = 1, limitj = iblockSize - 1;
                if (leftid == MPI_PROC_NULL) {
                    j = 2;
                }
                if (rightid == MPI_PROC_NULL) {
                    limitj = iblockSize - 2;
                }

                for (; j < limitj; j++) {
                    // top border cells
                    if (i > 1 || (upid != MPI_PROC_NULL && i == 1)) {
                        arrayB[i * iblockSize + j] = (arrayA[i * iblockSize + j] + arrayA[(i - 1) * iblockSize + j] +
                                                      arrayA[(i + 1) * iblockSize + j]
                                                      + arrayA[i * iblockSize + j - 1] +
                                                      arrayA[i * iblockSize + j + 1]) / 5;
                    }
                    // bottom border cells
                    if (i > 1 || (downid != MPI_PROC_NULL && i == 1)) {
                        arrayB[jblockSize * iblockSize - (i + 1) * iblockSize + j] =
                                (arrayA[jblockSize * iblockSize - (i + 1) * iblockSize + j]
                                 + arrayA[jblockSize * iblockSize - (i + 2) * iblockSize + j] +
                                 arrayA[jblockSize * iblockSize - (i + 1) * iblockSize + j + 1]
                                 + arrayA[jblockSize * iblockSize - (i) * iblockSize + j] +
                                 arrayA[jblockSize * iblockSize - (i + 1) * iblockSize + j - 1]) / 5;
                    }
                }
            }
            for (int j = 5; j < jblockSize - 5; j++) {
                for (int i = 1; i < 5; i++) {
                    // left border cells
                    if (i > 1 || (leftid != MPI_PROC_NULL && i == 1)) {
                        arrayB[j * iblockSize + i] = (arrayA[j * iblockSize + i] + arrayA[(j - 1) * iblockSize + i] +
                                                      arrayA[(j + 1) * iblockSize + i]
                                                      + arrayA[j * iblockSize + i - 1] +
                                                      arrayA[j * iblockSize + i + 1]) / 5;
                    }
                    // right border cells
                    if (i > 1 || (rightid != MPI_PROC_NULL && i == 1)) {
                        arrayB[(j + 1) * iblockSize - 1 - i] =
                                (arrayA[(j + 1) * iblockSize - 1 - i] + arrayA[(j) * iblockSize - 1 - i] +
                                 arrayA[(j + 2) * iblockSize - 1 - i]
                                 + arrayA[(j + 1) * iblockSize - 2 - i] + arrayA[(j + 1) * iblockSize - i]) / 5;
                    }
                }
            }
            // 2nd iteration
            for (int i = 2; i < 4; i++) {
                for (int j = 2; j < iblockSize - 2; j++) {
                    // top border cells
                    arrayA[i * iblockSize + j] = (arrayB[i * iblockSize + j] + arrayB[(i - 1) * iblockSize + j] +
                                                  arrayB[(i + 1) * iblockSize + j]
                                                  + arrayB[i * iblockSize + j - 1] + arrayB[i * iblockSize + j + 1]) /
                                                 5;
                    // bottom border cells
                    arrayA[jblockSize * iblockSize - (i + 1) * iblockSize + j] =
                            (arrayB[jblockSize * iblockSize - (i + 1) * iblockSize + j]
                             + arrayB[jblockSize * iblockSize - (i + 2) * iblockSize + j] +
                             arrayB[jblockSize * iblockSize - (i + 1) * iblockSize + j + 1]
                             + arrayB[jblockSize * iblockSize - (i) * iblockSize + j] +
                             arrayB[jblockSize * iblockSize - (i + 1) * iblockSize + j - 1]) / 5;
                    progress += std::abs(arrayB[i * iblockSize + j] - arrayA[i * iblockSize + j]);
                    progress += std::abs(arrayB[jblockSize * iblockSize - (i + 1) * iblockSize + j] -
                                         arrayA[jblockSize * iblockSize - (i + 1) * iblockSize + j]);
                }
            }
            for (int j = 4; j < jblockSize - 4; j++) {
                for (int i = 2; i < 4; i++) {
                    // left border cells
                    arrayA[j * iblockSize + i] = (arrayB[j * iblockSize + i] + arrayB[(j - 1) * iblockSize + i] +
                                                  arrayB[(j + 1) * iblockSize + i]
                                                  + arrayB[j * iblockSize + i - 1] + arrayB[j * iblockSize + i + 1]) /
                                                 5;
                    // right border cells
                    arrayA[(j + 1) * iblockSize - 1 - i] =
                            (arrayB[(j + 1) * iblockSize - 1 - i] + arrayB[(j) * iblockSize - 1 - i] +
                             arrayB[(j + 2) * iblockSize - 1 - i]
                             + arrayB[(j + 1) * iblockSize - 2 - i] + arrayB[(j + 1) * iblockSize - i]) / 5;
                    progress += std::abs(arrayB[j * iblockSize + i] - arrayA[j * iblockSize + i]);
                    progress += std::abs(arrayB[(j + 1) * iblockSize - 1 - i] - arrayA[(j + 1) * iblockSize - 1 - i]);
                }
            }
            // treat arrayB as send buffers
            if (upid != MPI_PROC_NULL)
                memcpy(&arrayB[0], &arrayA[2 * iblockSize], 2 * iblockSize * sizeof(double));
            if (downid != MPI_PROC_NULL)
                memcpy(&arrayB[iblockSize * jblockSize - 2 * iblockSize],
                       &arrayA[iblockSize * jblockSize - 4 * iblockSize], 2 * iblockSize * sizeof(double));

            for (int i = 0; i < jblockSize - 4; i++) {
                if (leftid != MPI_PROC_NULL)
                    memcpy(&arrayB[(i + 2) * iblockSize], &arrayA[(i + 2) * iblockSize + 2], 2 * sizeof(double));
                if (rightid != MPI_PROC_NULL)
                    memcpy(&arrayB[(i + 3) * iblockSize - 2], &arrayA[(i + 3) * iblockSize - 4], 2 * sizeof(double));
            }
            // send out/receive ghost-cells
            MPI_Request req[8];
            MPI_Isend(&arrayB[0], 2 * iblockSize, MPI_DOUBLE, upid, 0, MPI_COMM_WORLD, &req[0]); // up border
            MPI_Isend(&arrayB[iblockSize * jblockSize - 2 * iblockSize], 2 * iblockSize, MPI_DOUBLE, downid, 0,
                      MPI_COMM_WORLD, &req[1]);  //down border
            MPI_Isend(&arrayB[(3) * iblockSize - 2], 1, sideBorderType, rightid, 0, MPI_COMM_WORLD,
                      &req[2]);  //right border
            MPI_Isend(&arrayB[(2) * iblockSize], 1, sideBorderType, leftid, 0, MPI_COMM_WORLD, &req[3]);  //left border

            MPI_Irecv(&arrayA[(2) * iblockSize], 1, sideBorderType, leftid, 0, MPI_COMM_WORLD, &req[6]);
            MPI_Irecv(&arrayA[(3) * iblockSize - 2], 1, sideBorderType, rightid, 0, MPI_COMM_WORLD, &req[7]);
            MPI_Irecv(&arrayA[0], 2 * iblockSize, MPI_DOUBLE, upid, 0, MPI_COMM_WORLD, &req[4]);
            MPI_Irecv(&arrayA[iblockSize * jblockSize - 2 * iblockSize], 2 * iblockSize, MPI_DOUBLE, downid, 0,
                      MPI_COMM_WORLD, &req[5]);

            // calculate interiors
            for (int i = 5; i < jblockSize - 5; i++) {
                for (int j = 5; j < iblockSize - 5; j++) {
                    arrayB[i * iblockSize + j] = (arrayA[i * iblockSize + j] + arrayA[(i - 1) * iblockSize + j] +
                                                  arrayA[(i + 1) * iblockSize + j]
                                                  + arrayA[i * iblockSize + j - 1] + arrayA[i * iblockSize + j + 1]) /
                                                 5;
                }
            }
            for (int i = 4; i < jblockSize - 4; i++) {
                for (int j = 4; j < iblockSize - 4; j++) {
                    arrayA[i * iblockSize + j] = (arrayB[i * iblockSize + j] + arrayB[(i - 1) * iblockSize + j] +
                                                  arrayB[(i + 1) * iblockSize + j]
                                                  + arrayB[i * iblockSize + j - 1] + arrayB[i * iblockSize + j + 1]) /
                                                 5;
                    progress += std::abs(arrayA[i * iblockSize + j] - arrayB[i * iblockSize + j]);
                }
            }

            MPI_Allreduce(&progress, &reducedProgress, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

            if (myid == 0)
                iter += 2;

            // synchronize 
            MPI_Waitall(8, req, MPI_STATUSES_IGNORE);
        } while (reducedProgress >= epsilon);

        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::high_resolution_clock::now() - start);
        if (myid == 0)
            cout << elapsed.count() << endl;
//        if (output) {
//            int i = 0;
//            if(myid==0){
//                cout << "Iterations: " << iter << endl << endl;
//                print2Darray(arrayA, 0, 0, jblockSize, iblockSize, myid);
//                MPI_Send(&i,1,MPI_INTEGER,myid+1,0,MPI_COMM_WORLD);
//            }else if(myid != worldSize-1){
//                MPI_Recv(&i,1,MPI_INTEGER,myid-1,0,MPI_COMM_WORLD,MPI_STATUSES_IGNORE);
//                print2Darray(arrayA, 0, 0, jblockSize, iblockSize, myid);
//                MPI_Send(&i,1,MPI_INTEGER,myid+1,0,MPI_COMM_WORLD);
//            }else{
//                MPI_Recv(&i,1,MPI_INTEGER,myid-1,0,MPI_COMM_WORLD,MPI_STATUSES_IGNORE);
//                print2Darray(arrayA, 0, 0, jblockSize, iblockSize, myid);
//            }
//        }
        if (output)
            if (myid == 0)
                cout << "Iterations: " << iter << endl << endl;

        MPI_Type_free(&sideBorderType);
        delete[] arrayA;
        delete[] arrayB;
    }

    MPI_Finalize();

    return EXIT_SUCCESS;
}
