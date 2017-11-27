#include <iostream>
#include <cstdlib>


using namespace std;

class Wrapper {
public:
    float *arrayA;
    float *arrayB;
    int dim;
    int size;
    
    Wrapper(int dim, int size);
    ~Wrapper();
    void setBoundary1D(float left, float right);
    void setBoundary2D(float left, float right, float up, float down);
    void setBoundary3D(float left, float right, float up, float down, float front, float back);
    void print2DA();
    
};

Wrapper::Wrapper(int dim, int size) :size(size), dim(dim){
    this->dim = dim;
    this->size = size;
    
    //+2 to include boundary and therefore ignore border cases in jacobi
    if (dim == 1) {
        arrayA = new float[size + 2];
        arrayB = new float[size + 2];
    } else if (dim == 2) {
        arrayA = new float[(size +2) * (size +2)];
        arrayB = new float[(size +2) * (size +2)];
    } else if (dim == 3){
        arrayA = new float[(size +2) * (size +2) * (size +2)];
        arrayB = new float[(size +2) * (size +2) * (size +2)];
    } else {
        cerr << "wrong dimension" << endl;
    }
}

void Wrapper::setBoundary1D(float left, float right) {
    arrayA[0] = left;
    arrayA[size + 1] = right;
    
    arrayB[0] = left;
    arrayB[size + 1] = right;
}

void Wrapper::setBoundary2D(float left, float right, float up, float down) {
    for (int i = 0; i < size + 2; i++) {
        arrayA[i] = up;
        arrayA[i * (size + 2) + 0] = left;
        arrayA[i * (size + 2) + size + 1] = right;
        arrayA[(size + 1) * (size + 2) + i] = down;
        
        arrayB[i] = up;
        arrayB[i * (size + 2) + 0] = left;
        arrayB[i * (size + 2) + size + 1] = right;
        arrayB[(size + 1) * (size + 2) + i] = down;
    }
}

void Wrapper::setBoundary3D(float left, float right, float up, float down, float front, float back) {
    
}

void Wrapper::print2DA() {
    cout << "dim: " << this->dim << " size: " << this->size << endl;
    for (int i = 0; i < size + 2; i++) {
        for (int j = 0; i < size + 2; j++) {
            cout << arrayA[i * (size + 2) + j] << ", ";
        }
        cout << endl;
    }
}

Wrapper::~Wrapper() {
    delete arrayA;
    delete arrayB;
}

int main () {
    Wrapper test(2, 1);
    test.setBoundary2D(1,2,3,4);
    //cout << "dim: " << test.dim << " size: " << test.size << endl;
    test.print2DA();
}
