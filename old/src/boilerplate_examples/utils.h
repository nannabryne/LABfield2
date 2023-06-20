#include <stdlib.h>
#include "LATfield2.hpp"
using namespace LATfield2;


const string OUTPUT_PATH      = "../results/fresh_output/";
const string ORG_RESULTS_PATH = "../results/original_results/";



extern struct DiagnosticOutput{
    
    /* What kind of error? */
    bool error_is_relative; 
    bool error_is_absolute = error_is_relative ? false : true;

    /* Deviation from original sample data */
    double max_error;   // maximal error
    double min_error;   // minimal error
    double avg_error;   // average error

    /* Runtime values */
    double org_runtime;                 // original duration of computation in ms
    double new_runtime = org_runtime;   // duration of computation now in ms
    double enhancement = 100. * (org_runtime - new_runtime)
                        / org_runtime;  // performance enhancement

} diagnostics;


string file(const string filename, bool original=false){
    string file;
    if(original){
        COUT << endl << "Over-writing original results." << endl;
        file = ORG_RESULTS_PATH + filename;
    }
    else{
        file = OUTPUT_PATH + filename;
    }

    file += ".h5";

    return file;
}
