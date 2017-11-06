
#include <cstdlib>
#include <iostream>
#include <cstring>
#include <omp.h>
#include <ctime>

#include "chrono_timer.h"

using namespace std;

#define SMALL    32  // Arrays size <= SMALL switches to insertion sort

void merge(int a[], int size, int temp[]) {
    int i1 = 0;
    int i2 = size/2;
    int tempi = 0;
    while (i1 < size/2 && i2 < size) {
        if (a[i1] < a[i2]) {
            temp[tempi] = a[i1];
            i1++;
        } else {
            temp[tempi] = a[i2];
            i2++;
        }
        tempi++;
    }
    while (i1 < size/2) {
        temp[tempi] = a[i1];
        i1++;
        tempi++;
    }
    while (i2 < size) {
        temp[tempi] = a[i2];
        i2++;
        tempi++;
    }
    
    // Copy sorted temp array into main array, a
    memcpy(a, temp, size*sizeof(int));
}


void insertion_sort(int a[], int size) {
    int i;
    for (i=0; i < size; i++) {
        int j, v = a[i];
        for (j = i - 1; j >= 0; j--) {
            if (a[j] <= v)
                break;
            a[j + 1] = a[j];
        }
        a[j + 1] = v;
    }
}

void mergesort_serial(int a[], int size, int temp[]) {
    if (size < SMALL) {
        insertion_sort(a, size);
        return;
    }
    mergesort_serial(a, size/2, temp);
    mergesort_serial(a + size/2, size - size/2, temp);
    
    merge(a, size, temp);
}

void mergesort_parallel(int a[], int size, int temp[], int threads) {
    if ( threads == 1) {
        mergesort_serial(a, size, temp);
    } else if (threads > 1) {
#pragma omp parallel sections num_threads(2)
        {
#pragma omp section
            {
                mergesort_parallel(a, size/2, temp, threads/2);}
#pragma omp section
            {
                mergesort_parallel(a + size/2, size - size/2, temp + size/2, threads - threads/2);}
        }
        
        merge(a, size, temp);
    } else {
        printf("Error: %d threads\n", threads);
        return;
    }
}
    
int main(int argc, char* argv[]) {
    if ( argc != 3 ) {
        cout << "Usage: array-size number-of-threads" << endl;
        return EXIT_FAILURE;
    }
    
    int size = atoi(argv[1]); // array size
    int threads = atoi(argv[2]); // number of threads
    
    
    omp_set_nested(1);
    if (omp_get_nested() != 1) {
        cout << "Nested omp not supported" << endl;
    }
    
    int* a    = (int*) malloc(sizeof(int)*size);
    int* temp = (int*) malloc(sizeof(int)*size);
    if (a == NULL || temp == NULL) {
        cout << "Array allocation error" << endl;
        return EXIT_FAILURE;
    }
    
    // Random array initialization
    int i;
    srand(time(NULL));
    for(i=0; i<size; i++) {
        a[i] = rand() % size;
    }
    
    {
        ChronoTimer t("Sort");
        mergesort_parallel(a, size, temp, threads);
    }
    
    //TODO Unit test....
    // Result check
    for(i=1; i<size; i++) {
        if (!(a[i-1] <= a[i])) {
            cout << "Failed check!" << endl;
            return EXIT_FAILURE;
        }
    }
    cout << "check successful" << endl;
    
    return EXIT_SUCCESS;
}
