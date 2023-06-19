/**
 * @file field_manipulation.hpp
 * @author Nanna Bryne
 * 
 * @brief Boilerplate example: create two fields, initialise them and compute the linear combination
 * 
 * @details Uses functions ..... from LATfield2
*/

#include "LATfield2.hpp"
using namespace LATfield2;

#define str_filename "field_manipulation"


int run_example(int n, int m){

    parallel.initialize(n, m);

    const int dim       = 3;            // #{dimensions}
    int halo_size       = 1;            // size of halo
    int lat_size[dim]   = {40,40,40};   // lattice size



    /* Define lattice */

    Lattice lat(dim, lat_size, halo_size);

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
        A2(x) = sin(x.coord(0)*(x.coord(1)+M_PI)) + x.coord(2)*cos(x.coord(2));
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


    // const string filename = "field_manipulation";

    bool ORIGINAL = false;
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



int compare_results(const double tolerance=1e-10){

    string outfile = OUTPUT_PATH + str_filename + ".h5";
    string orgfile = ORG_RESULTS_PATH + str_filename + ".h5";


    Field<Real> B_fresh;
    // B_fresh.loadHDF5(outfile);


    return 0;
}


