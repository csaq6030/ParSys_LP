#define CATCH_CONFIG_MAIN
#include "../catch.hpp"

#include "n_queens_cpp.cpp"

TEST_CASE("Test nQueens problem", "rec nQueen") {

	uint32_t results[] = {1,0,0,2,10,4,40,92,352,724};

	for(uint32_t i = 1; i <=10; i++){
		setTest(i);
		CHECK(rec_n_queens() == results[i-1]);
	}
}

