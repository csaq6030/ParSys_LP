#ifndef WRAPPER_H
#define WRAPPER_H

#include <iostream>
#include <cstdlib>


using namespace std;

class Wrapper {
private:
    void print1DA();
    void print1DB();
    void print2DA();
    void print3DA();
    void print2DB();
    void print3DB();
public:
    float *arrayA;
    float *arrayB;
    const int dim;
    const int size;
    
    Wrapper(const int dim,const int size);
    ~Wrapper();
    void setBoundary(const float left, const float right);
    void setBoundary(const float left, const float right, const float up, const float down);
    void setBoundary(const float left, const float right, const float up, const float down, const float front, const float back);
    void printA();
    void printB();
    void print();
    
};

Wrapper::Wrapper(const int dim,const int size) :size(size), dim(dim){
    //+2 to include boundary and therefore ignore border cases in jacobi
    if (dim == 1) {
        arrayA = new float[size + 2]();
        arrayB = new float[size + 2]();
    } else if (dim == 2) {
        arrayA = new float[(size +2) * (size +2)]();
        arrayB = new float[(size +2) * (size +2)]();
    } else if (dim == 3){
        arrayA = new float[(size +2) * (size +2) * (size +2)]();
        arrayB = new float[(size +2) * (size +2) * (size +2)]();
    } else {
        cerr << "wrong dimension" << endl;
    }
}

//1D
void Wrapper::setBoundary(const float left, const float right) {
    arrayA[0] = left;
    arrayA[size + 1] = right;
    
    arrayB[0] = left;
    arrayB[size + 1] = right;
}

//2D
void Wrapper::setBoundary(const float left, const float right, const float up, const float down) {
    for (int i = 0; i < size + 2; i++) {
        arrayA[i] = up;
        arrayA[i * (size + 2)] = left;
        arrayA[i * (size + 2) + size + 1] = right;
        arrayA[(size + 1) * (size + 2) + i] = down;
        
		arrayB[i] = up;
        arrayB[i * (size + 2)] = left;
        arrayB[i * (size + 2) + size + 1] = right;
        arrayB[(size + 1) * (size + 2) + i] = down;
    }
}

//3D
void Wrapper::setBoundary(const float left, const float right, const float up, const float down, const float front, const float back) {
    for (int i = 0; i < size + 2; i++) {
        for (int j = 0; j < size + 2; j++) {
            arrayA[i * (size + 2) + j] = front;
            arrayA[(size + 1) * (size + 2) * (size + 2) + i * (size + 2) + j] = back;
            
            arrayB[i * (size + 2) + j] = front;
            arrayB[(size + 1) * (size + 2) * (size + 2) + i * (size + 2) + j] = back;
            }
    }
    
    //two loops because otherwise front would override left etc.
    for (int i = 0; i < size + 2; i++) {
        for (int j = 0; j < size + 2; j++) {
            arrayA[i *(size + 2) * (size + 2) + j * (size + 2)] = left;
            arrayA[i *(size + 2) * (size + 2) + j * (size + 2) + (size + 1)] = right;
            arrayA[i *(size + 2) * (size + 2) + j] = up;
            arrayA[i *(size + 2) * (size + 2) +(size + 1) * (size + 2) + j] = down;
            
            arrayB[i *(size + 2) * (size + 2) + j * (size + 2)] = left;
            arrayB[i *(size + 2) * (size + 2) + j * (size + 2) + (size + 1)] = right;
            arrayB[i *(size + 2) * (size + 2) + j] = up;
            arrayB[i *(size + 2) * (size + 2) +(size + 1) * (size + 2) + j] = down;
        }
    }
    
    
}

void Wrapper::print1DA() {
    cout << "Printing 1D A: " << endl;
    for (int i = 0; i < size + 2; i++) {
            cout << arrayA[i] << ", ";
        }
    cout << endl;
}

void Wrapper::print2DA() {
    cout << "Printing 2D A: " << endl;
    for (int i = 0; i < size + 2; i++) {
        for (int j = 0; j < size + 2; j++) {
            cout << arrayA[i * (size + 2) + j] << ", ";
        }
        cout << endl;
    }
}

void Wrapper::print3DA() {
    cout << "Printing 3D A: " << endl;
    for (int i = 0; i < size + 2; i++) {
        for (int j = 0; j < size + 2; j++) {
            for (int k = 0; k < size + 2; k++) {
               cout << arrayA[i * (size + 2) * (size + 2) + j * (size + 2) + k] << ", ";
            }
           cout << endl;
        }
        cout << endl;
    }
}

void Wrapper::print1DB() {
    cout << "Printing 1D B: " << endl;
    for (int i = 0; i < size + 2; i++) {
        cout << arrayB[i] << ", ";
    }
    cout << endl;
}

void Wrapper::print2DB() {
    cout << "Printing 2D B: " << endl;
    for (int i = 0; i < size + 2; i++) {
        for (int j = 0; j < size + 2; j++) {
            cout << arrayB[i * (size + 2) + j] << ", ";
        }
        cout << endl;
    }
}

void Wrapper::print3DB() {
    cout << "Printing 3D B: " << endl;
    for (int i = 0; i < size + 2; i++) {
        for (int j = 0; j < size + 2; j++) {
            for (int k = 0; k < size + 2; k++) {
                cout << arrayB[i * (size + 2) * (size + 2) + j * (size + 2) + k] << ", ";
            }
            cout << endl;
        }
        cout << endl;
    }
}

void Wrapper::printA() {
    if (dim == 1) {
        print1DA();
    } else if (dim == 2) {
        print2DA();
    } else {
        print3DA();
    }
}

void Wrapper::printB() {
    if (dim == 1) {
        print1DB();
    } else if (dim == 2) {
        print2DB();
    } else {
        print3DB();
    }
}

void Wrapper::print() {
    if (dim == 1) {
        print1DA();
        cout << endl << endl;
        print1DB();
    } else if (dim == 2) {
        print2DA();
        cout << endl << endl;
        print2DB();
    } else {
        print3DA();
        cout << endl << endl;
        print3DB();
    }
}

Wrapper::~Wrapper() {
    delete[] arrayA;
    delete[] arrayB;
}

#endif

/*
int main () {
    int size = 5;
    int dim = 3;
    Wrapper test(dim, size);
    test.setBoundary(1,2,3,4,5,6);
    test.arrayA[3 * (size + 2) * (size + 2) + 3 * (size + 2) + 3] = 42.42;
    test.printA();
} */
