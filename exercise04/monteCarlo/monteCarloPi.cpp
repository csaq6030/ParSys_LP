#include <iostream>
#include <random>
#include <cmath>
#include <chrono>
#include <omp.h>

#include "chrono_timer.h"

using namespace std;

double MCS_Seq(unsigned long n){
	unsigned long c = 0;
	double x, y;
	random_device rd;
	mt19937 gen(rd()); 
	uniform_real_distribution<double> dis(0, 1.0);
	for(unsigned i = 0; i < n; i++){
		x = dis(gen);
		y = dis(gen);
		if(sqrt(x*x + y*y) <= 1)
			c++;
	}
	return (double) c/n*4.0;
}

double MCS_Par(unsigned long n){
	unsigned long c = 0;
	#pragma omp parallel
	{
		double x, y;
		random_device rd;
		mt19937 gen(rd() + omp_get_thread_num()); 
		uniform_real_distribution<double> dis(0, 1.0);
		#pragma omp for reduction(+:c)
		for(unsigned i = 0; i < n; i++){
			x = dis(gen);
			y = dis(gen);
			if(sqrt(x*x + y*y) <= 1)
				c++;
		}	
	}
	return (double) c/n*4.0;
}

int main(int argc, char** argv) {

	if(argc != 3) return EXIT_FAILURE;
	unsigned long n = atoi(argv[1]);
	if(n == 0) return EXIT_FAILURE;
	string method(argv[2]); 

	double pi;
	
	if(!method.compare("seq")){
		//ChronoTimer("Pi Sequential");
		chrono::time_point<std::chrono::high_resolution_clock> start(std::chrono::high_resolution_clock::now());
    	pi = MCS_Seq(n);
		auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start);
		std::cout << "Timer - " << ": " << elapsed.count() << " ms " << std::endl;
    }else if(!method.compare("par")){
		chrono::time_point<std::chrono::high_resolution_clock> start(std::chrono::high_resolution_clock::now());
    	pi = MCS_Par(n);
		auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start);
		std::cout << "Timer - " << ": " << elapsed.count() << " ms " << std::endl;
	}else{
		cout << "executeable #OfSamples par/seq" << endl;
	}

	cout.precision(12);
	cout<< "Estimated value of pi is: " << pi << endl;
	cout<< "Actual    Value of pi is: 3.14159265359" << endl;

	return 0;
}









