#include <iostream>
#include <cstdlib>
#include <chrono>
#include <omp.h>

#include "solver.h"
#include "wrapper.h"

using namespace std;

void printError() {
    cout << "Usage: [number of threads] [dimension] [size] [epsilon] [left] [right] [up] [down] [front] [back]" << endl
    << "up and down are needed only for 2D and 3D" << endl
    << "front and back only for 3D" << endl;
}


int main(int argc, char* argv[]) {
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
		omp_set_num_threads(atoi(argv[1]));

        dim = atoi(argv[2]);
        size = atoi(argv[3]);
        epsilon = atof(argv[4]);
        left = atof(argv[5]);
        right = atof(argv[6]);
    }
    if (atoi(argv[2]) == 1 && argc == 7) { //1D
        Wrapper w1(1, size);
		w1.setBoundary(left, right);
		Solver s1(&w1, epsilon);
		
		const std::chrono::time_point<std::chrono::high_resolution_clock> start(std::chrono::high_resolution_clock::now());
		s1.compute();
		auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start);
		cout << elapsed.count() << endl;
        
    } else if (atoi(argv[2]) == 2 && argc == 9) { //2D
        const float up = atof(argv[7]);
        const float down = atof(argv[8]);

		Wrapper w1(2, size);
		w1.setBoundary(left, right, up, down);
		Solver s1(&w1, epsilon);
		
		const std::chrono::time_point<std::chrono::high_resolution_clock> start(std::chrono::high_resolution_clock::now());
		s1.compute();
		auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start);
		cout << elapsed.count() << endl;

		//s1.showResult();
        
    } else if (atoi(argv[2]) == 3 && argc == 11) { //3D
        const float up = atof(argv[7]);
        const float down = atof(argv[8]);
        const float front = atof(argv[9]);
        const float back = atof(argv[10]);

		Wrapper w1(3, size);
		w1.setBoundary(left, right, up, down, front, back);
		Solver s1(&w1, epsilon);
		
		const std::chrono::time_point<std::chrono::high_resolution_clock> start(std::chrono::high_resolution_clock::now());
		s1.compute();
		auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start);
		cout << elapsed.count() << endl;
        
        //cout <<"num_t: " << num_threads << " dim: " << dim << " size: " << size << " epsilon: " << epsilon << " left: " << left << " right: " << right << " up: " << up << " down: " << down << " front: " << front << " back: " << back << endl;
        
    } else {
        printError();
        return EXIT_FAILURE;
    }
}
