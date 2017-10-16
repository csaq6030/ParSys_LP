#include <vector>
#include <cstdlib>

using Matrix = std::vector<std::vector<double>>;

// initializes a square identity matrix of size n x n
Matrix id(unsigned n) {
	Matrix res;
	res.resize(n);

	//#pragma omp parallel for
	for(unsigned i=0; i<n; i++) {
		res[i].resize(n,0);
		res[i][i] = 1;
	}
	return res;
}

// computes the product of two matrices
Matrix operator*(const Matrix& a, const Matrix& b) {
	unsigned n = a.size();
	Matrix c = id(n);
	double tmp = 0;
	#pragma omp parallel for collapse(2) private(tmp) shared(a,b,c)
	for(unsigned i=0; i<n; ++i) {
		for(unsigned j=0; j<n; ++j) {
			tmp = 0;
			for(unsigned k=0; k<n; ++k) {
				tmp += a[i][k] * b[k][j];
			}
			c[i][j] = tmp;
		}
	}
	return c;
}


int main(int argc, char** argv) {
	
	if(argc!=2) return EXIT_FAILURE;
	unsigned n = atoi(argv[1]);
	if(n==0) return EXIT_FAILURE;

	// create two matrices
	auto a = id(n);
	a[0][0] = 42;
	auto b = id(n);

	// compute the product
	auto c = a * b;

	// check that the result is correct
	return (c == a) ? EXIT_SUCCESS : EXIT_FAILURE;
}
