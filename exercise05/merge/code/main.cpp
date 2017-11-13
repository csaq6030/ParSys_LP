#include "merge.cpp"

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
    srand(time(NULL));
    for(int i=0; i<size; i++) {
        a[i] = rand() % size;
    }
    
    {
        ChronoTimer t("Sort");
        mergesort_parallel(a, size, temp, threads);
    }
    
    return EXIT_SUCCESS;
}
