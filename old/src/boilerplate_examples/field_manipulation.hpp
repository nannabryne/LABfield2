/**
 * @file field_manipulation.hpp
 * @author Nanna Bryne
 * 
 * @brief Boilerplate example: 
 *          (1) create two fields A1 and A2 & initialise them
 *          (2) compute a linear combination B
 *          (3) save B
 *          (4) compare the new results to original results
 * 
 * @details Uses functions ..... from LATfield2
*/


#include "utils.h"


#define str_filename "field_manipulation"


struct DiagnosticOutput diagnostics;

struct Parameters{
    const int dim           = 3;            // #{dimensions}
    const int halo_size     = 1;            // size of halo
    const int lat_size[3]   = {40,40,40};   // lattice size
} params;


int run_example(int n, int m){

    parallel.initialize(n, m);

    /* Define lattice */

    Lattice lat(params.dim, params.lat_size, params.halo_size);

    /* Set up fields */
    Field<Real> A1;
    A1.initialize(lat);
    A1.alloc();
    Field<Real> A2;
    A2.initialize(lat);
    A2.alloc();

    Site x(lat);

    for(x.first(); x.test(); x.next()){
        A1(x) = sin(x.coord(0) + x.coord(2)/M_PI) + cos(x.coord(1));
        A2(x) = sin(x.coord(0)*(x.coord(1)+1.*M_PI)) + x.coord(2)*cos(x.coord(2));
    }

    A1.updateHalo();
    A2.updateHalo();

    /* Get a linear combination of the two fields */

    double alpha1 = 2.;
    double alpha2 = 1.2;

    Field<Real> B;
    B.initialize(lat);
    B.alloc();

    for(x.first(); x.test(); x.next()){
        B(x) = alpha1*A1(x) + alpha2*A2(x);
    }

    B.updateHalo();

    parallel.barrier();


    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


    const bool ORIGINAL = false;
    string file;
    if(ORIGINAL){
        COUT << endl << "Over-writing original results." << endl;
        file = ORG_RESULTS_PATH + str_filename;
    }
    else{
        file = OUTPUT_PATH + str_filename;
    }


#ifdef HDF5
    file += ".h5";
    B.saveHDF5(file);
#else
    file += ".txt";
    B.write(file);
#endif


    return 0;
}



void compare_results(const double tolerance=1e-10){

    /* Read files and create fields */

    Lattice lat(params.dim, params.lat_size, params.halo_size);

    string outfile = OUTPUT_PATH + str_filename + ".h5";
    string orgfile = ORG_RESULTS_PATH + str_filename + ".h5";


    Field<Real> B_fresh;
    B_fresh.initialize(lat); 
    B_fresh.alloc();
    B_fresh.loadHDF5(outfile);

    Field<Real> B_org;
    B_org.initialize(lat); 
    B_org.alloc();
    B_org.loadHDF5(orgfile);

    /* Compare fields 
    (going with absolute error instead of relative error to avoid division by zero) */

    Site x(lat);

    double abs_err;         // absolute error
    double max_err = 0.;    // maximal error
    double min_err = 20.;   // minimal error
    double avg_err = 0.;    // average error

    for(x.first(); x.test(); x.next()){
        abs_err = fabs( B_fresh(x) - B_org(x) );
        avg_err += abs_err;
        if(min_err > abs_err) min_err = abs_err;
        if(max_err < abs_err) max_err = abs_err;   
    }

    parallel.max(max_err);
    parallel.min(min_err);
    parallel.sum(avg_err);

    avg_err /= (double)lat.sites();

    parallel.barrier();
    
    COUT << max_err << endl;

    
    diagnostics.error_is_relative = true;

    diagnostics.max_error = max_err;
    diagnostics.min_error = min_err;
    diagnostics.avg_error = avg_err;
    diagnostics.org_runtime = 19.5;
    
}


