#ifndef HELPER_H
#define HELPER_H

#include <cmath>
#include <iostream>

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

#endif