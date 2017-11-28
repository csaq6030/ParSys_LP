#ifndef SOLVER_H
#define SOLVER_H


#include "wrapper.h"

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>

using namespace std;

class Solver {
public:
	Wrapper *data;
	const float epsilon;
	size_t iter;
	const int size;   

	vector<float*> indirectionA;
	vector<float*> indirectionB;

    Solver(Wrapper* d, const float epsilon);
	void compute();
	void showResult();
private:
	float iteration1D_A();	//input A-array, output B-array
	float iteration1D_B();
	float iteration2D_A();
	float iteration2D_B();
	float iteration3D_A();
	float iteration3D_B();
	int access(int i, int j);
	int access(int i, int j, int k);
    
};

Solver::Solver(Wrapper * d, const float epsilon) :data(d), epsilon(epsilon), size(data->size){
	if(data->dim == 2){
		for(unsigned i = 0; i < size+2; ++i) {
		indirectionA.resize(size+2);
		indirectionB.resize(size+2);
			for(unsigned j = 0; j < size+2; ++j) {
				indirectionA[i] = &data->arrayA[i*(size+2)];
				indirectionB[i] = &data->arrayB[i*(size+2)];
			}
		}
	}
}

void Solver::compute(){
	float progress = 500000;
	iter=0;
	switch(data->dim){
		case 1:
			while(progress > epsilon){
				if(iter%2 == 0){
					progress = iteration1D_A();
				} else{
					progress = iteration1D_B();
				}
				iter++;
			}			
			break;
		case 2:
			while(progress > epsilon){
				if(iter%2 == 0){
					progress = iteration2D_A();
				} else{
					progress = iteration2D_B();
				}
				iter++;
			}	
			break;
		case 3:
			while(progress > epsilon){
				if(iter%2 == 0){
					progress = iteration3D_A();
				} else{
					progress = iteration3D_B();
				}
				iter++;
			}	
			break;
		default:
			break;
	}
}

void Solver::showResult(){
	if(iter%2 == 0)
		data->printA();
	else
		data->printB();
	cout << "Iterations = " << iter << endl;
}

float Solver::iteration1D_A(){
	float prog = 0;
	#pragma omp parallel for simd reduction(+: prog)
	for(int i = 1; i < size + 1; i++){
		float tmp = (data->arrayA[i-1] + data->arrayA[i] + data->arrayA[i+1])/3;
		prog += abs(data->arrayA[i]-tmp);
		data->arrayB[i] = tmp;
	}
	return prog;
}

float Solver::iteration1D_B(){
	float prog = 0;
	#pragma omp parallel for simd reduction(+: prog)
	for(int i = 1; i < size + 1; i++){
		float tmp = (data->arrayB[i-1] + data->arrayB[i] + data->arrayB[i+1])/3;
		prog += abs(data->arrayB[i]-tmp);
		data->arrayA[i] = tmp;
	}
	return prog;
}

float Solver::iteration2D_A(){
	float prog = 0;
	#pragma omp parallel for simd reduction(+: prog)
	for(int i = 1; i < size + 1; i++){
		for(int j = 1; j < size + 1; j++){
			float tmp = (indirectionA[i][j-1] + indirectionA[i][j] + indirectionA[i][j+1] + indirectionA[i-1][j] + indirectionA[i+1][j])/5;
			prog += abs(indirectionA[i][j]-tmp);
			data->arrayB[access(i, j)] = tmp;
		}
	}
	return prog;
}

float Solver::iteration2D_B(){
	float prog = 0;
	#pragma omp parallel for simd reduction(+: prog)
	for(int i = 1; i < size + 1; i++){
		for(int j = 1; j < size + 1; j++){
			float tmp = (indirectionB[i][j-1] + indirectionB[i][j] + indirectionB[i][j+1] + indirectionB[i-1][j] + indirectionB[i+1][j])/5;
			prog += abs(indirectionB[i][j]-tmp);
			data->arrayA[access(i, j)] = tmp;
		}
	}
	return prog;
}

float Solver::iteration3D_A(){
	float prog = 0;
	#pragma omp parallel for reduction(+: prog) 
	for(int i = 1; i < size + 1; i++){
		for(int j = 1; j < size + 1; j++){
			for(int k = 1; k < size + 1; k++){
				float tmp = (data->arrayB[access(i-1, j, k)] + data->arrayB[access(i+1, j, k)] + data->arrayB[access(i, j, k)] +
							 data->arrayB[access(i, j-1, k)] + data->arrayB[access(i, j+1, k)] + data->arrayB[access(i, j, k-1)] + data->arrayB[access(i, j, k+1)] )/7;
				prog += abs(data->arrayB[access(i, j, k)]-tmp);
				data->arrayA[access(i, j, k)] = tmp;
			}			
		}
	}
	return prog;
}

float Solver::iteration3D_B(){
	float prog = 0;
	#pragma omp parallel for reduction(+: prog) 
	for(int i = 1; i < size + 1; i++){
		for(int j = 1; j < size + 1; j++){
			for(int k = 1; k < size + 1; k++){
				float tmp = (data->arrayA[access(i-1, j, k)] + data->arrayA[access(i+1, j, k)] + data->arrayA[access(i, j, k)] +
							 data->arrayA[access(i, j-1, k)] + data->arrayA[access(i, j+1, k)] + data->arrayA[access(i, j, k-1)] + data->arrayA[access(i, j, k+1)] )/7;
				prog += abs(data->arrayA[access(i, j, k)]-tmp);
				data->arrayB[access(i, j, k)] = tmp;
			}			
		}
	}
	return prog;
}

inline int Solver::access(int i, int j){
	return ((i) * (size + 2) + j);
}

inline int Solver::access(int i, int j, int k){
	return (i *(size + 2) * (size + 2) + j* (size + 2) + k);
}

#endif

/*
int main(){
	int size = 5;
    int dim = 3;
	float ep = 10;
    Wrapper test(dim, size);
    test.setBoundary(1,2,3,4,5,6);

	Solver s(&test,ep);
	s.compute();
    s.showResult();

	
}*/

