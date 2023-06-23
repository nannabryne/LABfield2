/**
 * @file fft_execution/main.cpp
 * @author Nanna Bryne
 * 
 * @brief Unit test: for FFT, both forward and backward
 * 
 * @details 
 *      - Only consider quadratic lattices
 * 
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
    const int dim           = 3;            // #{dimensions}
    const int halo          = 1;            // size of halo
    const int halo_k        = 0;            // size of halo in Fourier 'space'
    const int npts          = 512;         // lattice size in one dimension


    Diagnostics d("fft_execution", n, m);
    d.provide_reference_parameters(64, 2400.);
    d.provide_computation_parameters(npts);

    parallel.initialize(n, m);
    num_prcs = n*m;

    COUT << " " << endl;


    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


    /* Initialise lattices and sites */

    Lattice lat_x(dim, npts, halo);
    Lattice lat_k;
    lat_k.initializeRealFFT(lat_x, halo_k);

    Site x(lat_x);
    rKSite k(lat_k);

    /* Create fields */


    Field<Real> f_x; f_x.initialize(lat_x, 1);    // f(x)
    Field<Imag> f_k; f_k.initialize(lat_k, 1);    // f(k)

    PlanFFT<Imag> plan_f(&f_x, &f_k);   // FFT planner for f(x) -> f(k)


    double sigma_sqrd   = lat_x.size(0)*lat_x.size(1)/10; 
    double mu[3]        = {lat_x.size(0)*0.5, lat_x.size(1)*0.5, lat_x.size(2)*0.5};
    double A            = pow(2.*M_PI*sigma_sqrd, -dim*0.5); 

    double xp1, xp2, xp3;   // x' = x - μ = (x'_1, x'_2, x'_3) 

    for(x.first(); x.test(); x.next()){
        xp1 = (x.coord(0) - mu[0]); xp1 *= xp1;
        xp2 = (x.coord(1) - mu[1]); xp2 *= xp2;
        xp3 = (x.coord(2) - mu[2]); xp3 *= xp3;

        f_x(x) = A*exp(-(xp1+xp2+xp3) / sigma_sqrd*.5);
    }

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

    double timer_pause, timer_play;

    // ! START TIMING HERE
    timer_start = MPI_Wtime();

    COUT << endl << "FFT forward ..." << endl;
    timer_play = MPI_Wtime();
    plan_f.execute(FFT_FORWARD);
    timer_pause = MPI_Wtime();
    COUT << "finished in " << (timer_pause-timer_play)*1000.0 << " ms" << endl;


    COUT << endl << "FFT backward ..." << endl;
    timer_play = MPI_Wtime();
    plan_f.execute(FFT_BACKWARD);
    timer_pause = MPI_Wtime();
    COUT << "finished in " << (timer_pause-timer_play)*1000.0 << " ms" << endl;

    // ! END TIMING HERE
    runtime = ((MPI_Wtime() - timer_start)*1000.0);

    f_x.updateHalo();


    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


#ifdef _BENCH
    f_x.saveHDF5(orgfile);
    COUT << "reference runtime: " << runtime << " ms" << endl;
    COUT << "used " << num_prcs << " processes" << endl;
#else

    f_x.saveHDF5(outfile);

    
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

