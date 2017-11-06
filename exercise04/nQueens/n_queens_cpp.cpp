#include <cstdint>
#include <cstdlib>
#include <iostream>


#include "chrono_timer.h"

using namespace std;

uint32_t solution = 0;
uint32_t* ptr_s = &solution;
uint32_t N;

void putQueen(uint32_t *queens, uint32_t row, uint32_t col) {
    for (uint32_t i = 0; i < row; i++) {
        // check same row and diagonal
        if (queens[i] == col || queens[i] == col + row - i || queens[i] == col - row +i) {
            return;
        }
    }
    
    //put Queen
    queens[row] = col;
    
    //found a solution
    if(row == N-1) {
        //solution++;
        __sync_add_and_fetch(ptr_s, 1);
    } else {
        for(uint32_t i = 0; i < N; i++) {
            putQueen(queens, row+1, i);
        }
    }
}

void rec_n_queens() {
#pragma omp parallel for
    for(uint32_t i = 0; i < N; i++) {
        uint32_t *queens = (uint32_t*)calloc(N, sizeof(uint32_t));
        putQueen(queens, 0, i);
        free(queens);
    }
}

int main(int argc, char *argv[]) {
	if(argc == 2) {
		N = atoi(argv[1]);	
	} else {
        cout << "Usage: n_queens [N]" << endl;
		return EXIT_FAILURE;
	}	
    
    {
        ChronoTimer t("Execution ");

        rec_n_queens();
        
        
    
    }
    
    cout << "Found " << solution <<  " solutions" << endl ;

	return EXIT_SUCCESS;

}
