/**
 * @file field_manipulation/main.cpp
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


#include "../utils.h"


// struct DiagnosticOutput diagnostics;

struct Parameters{
    const int dim           = 3;            // #{dimensions}
    const int halo_size     = 1;            // size of halo
    const int lat_size[3]   = {40,40,40};   // lattice size
    const int N_lat         = 40;           // lattice size in one dimension
} params;




int main(int argc, char **argv){

    double time_ref, runtime;

    string orgfile = "org_output.h5";
    string outfile = "fresh_output.h5";

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

    Diagnostics d("field_manipulation", n, m);
    parallel.initialize(n, m);
    
    COUT << " " << endl;

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

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

    // ! START TIMING HERE
    time_ref = MPI_Wtime();

    for(x.first(); x.test(); x.next()){
        B(x) = alpha1*A1(x) + alpha2*A2(x);
    }

    // ! END TIMING HERE
    runtime = 1000.0 * (MPI_Wtime() - time_ref);

    B.updateHalo();

    parallel.barrier();


    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

    
    B.saveHDF5(outfile);

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

    /* Ensure validity of output output */

    parallel.barrier();
    field_comparison(params.dim, params.N_lat, params.halo_size, &d);

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

    /* Diagnostic analysis*/

    d.performance_metrics(0.5, runtime);

    d.write_epicrisis();
    d.print_epicrisis();


    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */



    COUT << " " << endl;


    return 0;
}
