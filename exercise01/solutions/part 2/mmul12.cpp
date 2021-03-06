#include <vector>
#include <cstdlib>

using Matrix = std::vector<std::vector<double>>;

// initializes a square identity matrix of size n x n
Matrix id(const unsigned n) {
	Matrix res = std::vector<std::vector<double>>(n, std::vector<double>(n));
	for(unsigned i=0; i<n; i++) {
		//resize by filling 0s to the vector, the position of the 1s in an identity already know -> remove 2nd iteration
		//res[i].resize(n,0);
		res[i][i] = 1;
		/*
		for(unsigned j=0; j<n; j++) {
			res[i][j] = (i == j) ? 1 : 0;
		}*/
	}
	return res;
}

// computes the product of two matrices
Matrix operator*(const Matrix& a, const Matrix& b) {
	const unsigned n = a.size();
    // no id Matrix needed
	Matrix c = std::vector<std::vector<double>>(n, std::vector<double>(n));
	// store result in tmp to reduce operator[] calls
	double tmp = 0;
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
	const unsigned n = atoi(argv[1]);
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
