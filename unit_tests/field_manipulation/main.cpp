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



int main(int argc, char **argv){

    double timer_start, runtime;

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

    // simulation parameters:
    const int dim   = 3;        // #{dimensions}
    const int halo  = 1;        // size of halo
    const int npts  = 512;     // lattice size in one dimension

    Diagnostics d("field_manipulation", n, m);
    d.provide_reference_parameters(64, 86.5);
    d.provide_computation_parameters(npts);
    


    parallel.initialize(n, m);
    
    COUT << " " << endl;

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

    /* Define lattice */

    Lattice lat(dim, npts, halo);

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
    timer_start = MPI_Wtime();

    for(x.first(); x.test(); x.next()){
        B(x) = alpha1*A1(x) + alpha2*A2(x);
    }

    // ! END TIMING HERE
    runtime = 1000.0 * (MPI_Wtime() - timer_start);

    B.updateHalo();

    parallel.barrier();


    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

    
#ifdef _BENCH
    B.saveHDF5(orgfile);
    COUT << "reference runtime: " << runtime << " ms" << endl;
    COUT << "used " << n*m << " processes" << endl;
#else
    
    B.saveHDF5(outfile);

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

    /* Ensure validity of output output */

    parallel.barrier();
    field_comparison(dim, npts, halo, &d);

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

    /* Diagnostic analysis*/


    d.compute_performance_metrics(runtime);

    d.write_epicrisis();
    d.print_epicrisis();

#endif
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
    


    COUT << " " << endl;


    return 0;
}
