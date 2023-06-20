/**
 * @file fft_execution.hpp
 * @author Nanna Bryne
 * 
 * @brief Boilerplate example: execute fast Fourier transform of XXX function
*/

#include "utils.h"

#define str_filename "fft_execution"

struct DiagnosticOutput diagnostics;

struct Parameters{
    const int dim           = 3;            // #{dimensions}
    const int halo_size     = 1;            // size of halo
    const int lat_size[3]   = {60,60,60};   // lattice size
    const int num_comp      = 1;            // #{components} ( = 1 for scalar field)

    const int halo_size_k   = 0;
} params;


int run_example(double *runtime){

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


    double time_ref, time_forward, time_backward, time_end;

    time_ref = MPI_Wtime();
    plan_f.execute(FFT_FORWARD);
    plan_f.execute(FFT_BACKWARD);
    *runtime = ((MPI_Wtime() - time_ref)*1000.0);


    f_x.updateHalo();


    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

    const bool ORIGINAL = false;
    
    f_x.saveHDF5(file(str_filename));
    

    return 0;
}
