#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "time_ms.h"

int s = 0;
int* ptr_s = &s;

static int is_safe(const int rows[],const int x,const int y) {
    
    if (y == 0)
        return 1;
    
    
    
    //much more faster then omp
    for (int i = 0; i < y; ++i) {  // same row and diagonals check
        if (rows[i] == x || rows[i] == x + y - i || rows[i] == x - y +i)
            return 0;
    }
    
    return 1;
    
    /*
    int out = 1;
    #pragma omp parallel for
    for (int i = 0; i < y; ++i) {  // same row and diagonals check
        if (rows[i] == x || rows[i] == x + y - i || rows[i] == x - y +i) {
            out = 0;
        }
    }
    return out; */
}

//unsafe with omp
static void draw_board(const int rows[],const int n) {
    printf("\nSolution #%d:\n---------------------------------\n", s);
    for (int y = 0; y < n; ++y) {
        for (int x = 0; x < n; ++x)
            printf(x == rows[y] ? "| Q " : "|   ");
        printf("|\n---------------------------------\n");
    }
}

void queens_helper(int rows[], int y,const int n) {
    int rows_local[n];
    #pragma omp parallel for private(rows_local)
    for (int x = 0; x < n; ++x) {
        memcpy(rows_local, rows, n*sizeof(int));
        if (is_safe(rows_local, x, y)) {
            rows_local[y] = x;
            if (y == n - 1) {
                __sync_add_and_fetch(ptr_s, 1);
                //draw_board(rows, n);
            }
            else
                queens_helper(rows_local, y+1, n);
        }
    }
}

int main(int argc, char** argv) {
    
    if(argc != 2) {
        printf("Usage: n_queens [N]\n");
        return EXIT_FAILURE;
    }
    
    const int n = atoi(argv[1]);
    int rows[n];
    
    unsigned long start_time = time_ms();
    queens_helper(rows, 0, n);
    printf("Found %i solutions\n", s);
    printf("Time: %9lu ms\n", time_ms() - start_time);
    
    return EXIT_SUCCESS;
}
