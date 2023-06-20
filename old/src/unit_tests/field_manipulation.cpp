/**
 * @file file_manipulation.hpp
 * @author Nanna Bryne
 * 
 * @brief Unit test: run code from '../boilerplate_examples/field_manipulation.hpp' and get diagnostic data
 * 
*/


#include "../boilerplate_examples/field_manipulation.hpp"


int main(int argc, char **argv){

    int n, m, num_prcs;

    for(int i=1; i<argc; i++){
        if(argv[i][0] != '-')
            continue;
        switch(argv[i][1]){
            case 'n':
                n = atoi(argv[++i]);
                break;
            case 'm':
                m = atoi(argv[++i]);
                break;
        }
    }
  

    // Start timing!
    double start = MPI_Wtime();

    run_example(n, m);
    
    // End timing!
    double stop = MPI_Wtime();

    COUT << " " << endl;

    /* Compare results */

    compare_results();




    /* Run diagnostics */

    double duration = (stop-start)*1000.0;

    COUT << "runtime: "<< duration << " ms"<< endl;
    diagnostics.new_runtime = duration;






    COUT << " " << endl;

    return 0;
}