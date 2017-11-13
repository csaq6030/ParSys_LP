#include <cstdint>
#include <cstdlib>
#include <iostream>


#include "chrono_timer.h"

using namespace std;

uint32_t solution = 0;
uint32_t* ptr_s = &solution;
uint32_t N;

void putQueen(uint32_t *queens, const uint32_t row, const uint32_t col) {
    for (uint32_t i = 0; i < row; i++) {
        // check same row and diagonal
        if (queens[i] == col || queens[i] == col + row - i || queens[i] == col - row +i) {
            return;
        }
    }
    
    //found a solution
    if(row == N-1) {
        //#pragma omp atomic
        //solution++;
        __sync_add_and_fetch(ptr_s, 1);
    } else {
        //put Queen
        queens[row] = col;
        
        
        for(uint32_t i = 0; i < N; i++) {
            putQueen(queens, row+1, i);
        }
    }
}

uint32_t rec_n_queens() {
    #pragma omp parallel
    {
        #pragma omp single nowait
        {
            for(uint32_t i = 0; i < N; i++) {
                #pragma omp task
                {
                    uint32_t queens[N];
                    putQueen(queens, 0, i);
                }
            }
        }
        #pragma omp taskwait
    }
	return solution;
}


void setTest(uint32_t n){
	N = n;
	solution = 0;
}
