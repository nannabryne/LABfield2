/**
 * @file fft_execution/main.cpp
 * @author Nanna Bryne
 * 
 * @brief Unit test: ...
 * 
*/


#include "../utils.h"


struct Parameters{
    const int dim           = 3;            // #{dimensions}
    const int halo_size     = 1;            // size of halo
    const int N_lat         = 60;           // lattice size in one dimension
    const int lat_size[3]   = {60,60,60};   // lattice size
    const int num_comp      = 1;            // #{components} ( = 1 for scalar field)

    const int halo_size_k   = 0;
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

    Diagnostics d("fft_execution", n, m);

    parallel.initialize(n, m);

    COUT << " " << endl;


    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


    /* Initialise lattices and sites */

    Lattice lat_x(params.dim, params.lat_size, params.halo_size);
    Lattice lat_k;
    lat_k.initializeRealFFT(lat_x, params.halo_size_k);

    Site x(lat_x);
    rKSite k(lat_k);

    /* Create fields */


    Field<Real> f_x; f_x.initialize(lat_x, params.num_comp);    // f(x)
    Field<Imag> f_k; f_k.initialize(lat_k, params.num_comp);    // f(k)

    PlanFFT<Imag> plan_f(&f_x, &f_k);   // FFT planner for f(x) -> f(k)


    double sigma_sqrd   = lat_x.size(0)*lat_x.size(1)/10; 
    double mu[3]        = {lat_x.size(0)*0.5, lat_x.size(1)*0.5, lat_x.size(2)*0.5};
    double A            = pow(2.*M_PI*sigma_sqrd, -params.dim*0.5); 

    double xp1, xp2, xp3;   // x' = x - μ = (x'_1, x'_2, x'_3) 

    for(x.first(); x.test(); x.next()){
        xp1 = (x.coord(0) - mu[0]); xp1 *= xp1;
        xp2 = (x.coord(1) - mu[1]); xp2 *= xp2;
        xp3 = (x.coord(2) - mu[2]); xp3 *= xp3;

        f_x(x) = A*exp(-(xp1+xp2+xp3) / sigma_sqrd*.5);
    }

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

    // ! START TIMING HERE
    time_ref = MPI_Wtime();

    plan_f.execute(FFT_FORWARD);
    plan_f.execute(FFT_BACKWARD);

    // ! END TIMING HERE
    runtime = ((MPI_Wtime() - time_ref)*1000.0);

    f_x.updateHalo();

    f_x.saveHDF5(outfile);


    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

    /* Ensure validity of output output */

    parallel.barrier();
    field_comparison(params.dim, params.N_lat, params.halo_size, &d);

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

    /* Diagnostic analysis*/

    d.performance_metrics(4.5, runtime);

    d.write_epicrisis();
    d.print_epicrisis();

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


    COUT << " " << endl;

    return 0;
}

