#include <iostream>
#include <cstdlib>

using namespace std;

void printError() {
    cout << "Usage: [number of threads] [dimension] [size] [epsilon] [left] [right] [up] [down] [front] [back]" << endl
    << "up and down are needed only for 2D and 3D" << endl
    << "front and back only for 3D" << endl;
}


int main(int argc, char* argv[]) {
    int num_threads = 0;
    int dim = 0;
    int size = 0;
    float epsilon = 0.;
    float left = 0;
    float right = 0;
    
    
    if ( argc > 11  || argc < 7) {
        printError();
        return EXIT_FAILURE;
    }
    
    if (argc >= 7) {
        num_threads = atoi(argv[1]);
        dim = atoi(argv[2]);
        size = atoi(argv[3]);
        epsilon = atof(argv[4]);
        left = atof(argv[5]);
        right = atof(argv[6]);
    }
    if (atoi(argv[2]) == 1 && argc == 7) { //1D
        
        
    } else if (atoi(argv[2]) == 2 && argc == 9) { //2D
        const float up = atof(argv[7]);
        const float down = atof(argv[8]);
        
    } else if (atoi(argv[2]) == 3 && argc == 11) { //3D
        const float up = atof(argv[7]);
        const float down = atof(argv[8]);
        const float front = atof(argv[9]);
        const float back = atof(argv[10]);
        
        cout <<"num_t: " << num_threads << " dim: " << dim << " size: " << size << " epsilon: " << epsilon << " left: " << left << " right: " << right << " up: " << up << " down: " << down << " front: " << front << " back: " << back << endl;
        
    } else {
        printError();
        return EXIT_FAILURE;
    }
}
