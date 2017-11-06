#define CATCH_CONFIG_MAIN
#include "../catch.hpp"

#include "merge.cpp"

TEST_CASE ("test merge", "[mergesort_parallel]") {
    int size = 1000;
    int threads = 4;
    
    int* a    = (int*) malloc(sizeof(int)*size);
    int* temp = (int*) malloc(sizeof(int)*size);
    
    // Random array initialization
    srand(time(NULL));
    for(int i=0; i<size; i++) {
        a[i] = rand() % size;
    }
    
    mergesort_parallel(a, size, temp, threads);
    
    CHECK (is_sorted(a, size));
}
