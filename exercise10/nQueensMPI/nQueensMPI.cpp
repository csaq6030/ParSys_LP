#include <cstdint>
#include <cstdlib>
#include <iostream>
#include "mpi.h"

#define TYPE uint32_t

using namespace std;

TYPE solution = 0;
TYPE N;


void putQueen(TYPE *queens, const TYPE row, const TYPE col) {
    
    for (TYPE i = 0; i < row; i++) {
        // check same row and diagonal
        if (queens[i] == col || queens[i] == col + row - i || queens[i] == col - row +i) {
            return;
        }
    }
    
    //found a solution
    if(row == N-1) {
        solution++;
    } else {
        //put Queen
        queens[row] = col;
        
        
        for(TYPE i = 0; i < N; i++) {
            putQueen(queens, row+1, i);
        }
    }
}

void rec_n_queens(TYPE start, TYPE end) {
    for(TYPE i = start; i < end; i++) {
        TYPE queens[N];
        putQueen(queens, 0, i);
    }
}

int main(int argc, char *argv[]) {
    int numtasks,taskid;
    double start, end;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);

    
    if(argc == 2) {
        N = atoi(argv[1]);
    } else {
        cout << "Usage: n_queens [N]" << endl;
        return EXIT_FAILURE;
    }
    
    start = MPI_Wtime();
    
    TYPE split = N / numtasks;
    
    if (taskid != numtasks - 1) { //needed for not devisable problem sizes
        rec_n_queens(taskid * split, taskid * split + split);
    } else {
        rec_n_queens(taskid * split, N);
    }
    
    TYPE reducedCount;
    MPI_Reduce(&solution, &reducedCount, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    
    end = MPI_Wtime();
    
    if(taskid == 0) {
        cout << end - start << endl;
        cout << "Results: " << reducedCount << endl;
    }
    
    
    MPI_Finalize();
    
    return EXIT_SUCCESS;
    
}


