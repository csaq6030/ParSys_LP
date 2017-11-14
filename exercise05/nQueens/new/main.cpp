#include "n_queens_cpp.cpp"

int main(int argc, char *argv[]) {
	if(argc == 2) {
		N = atoi(argv[1]);	
	} else {
        cout << "Usage: n_queens [N]" << endl;
		return EXIT_FAILURE;
	}	
    
    uint32_t result;
    {
        ChronoTimer t("Execution ");

         result = rec_n_queens();
    
    }
    
    cout << "Solutions;" << result <<  endl ;

	return EXIT_SUCCESS;

}
