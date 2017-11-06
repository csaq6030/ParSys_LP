#include "n_queens_cpp.cpp"

int main(int argc, char *argv[]) {
	if(argc == 2) {
		N = atoi(argv[1]);	
	} else {
        cout << "Usage: n_queens [N]" << endl;
		return EXIT_FAILURE;
	}	
    
    {
        ChronoTimer t("Execution ");

        rec_n_queens();
    
    }
    
    cout << "Found " << solution <<  " solutions" << endl ;

	return EXIT_SUCCESS;

}
