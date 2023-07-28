

#include "utils.h"

#ifdef LOOPCORR
    #include "test_LoopCorrection.hpp"
#endif
#ifdef PROJECTION
    #include "test_PartMeshProjection.hpp"
#endif
#ifdef FOURIER
    #include "test_FasterFourierTransform.hpp"
#endif
#ifdef PART
    #include "test_ParticleUpdate.hpp"
#endif


inline void initialiseParticles();

int main(int argc, char **argv){


    ///////////////////////////////////////////////////////////

    int n, m, o;  // parallel grid
#ifndef _OPENMP
    o = 1;
#else
    #pragma omp parallel
    {
        #pragma omp single
        o = omp_get_num_threads();
    }
#endif


    for(int i=1; i<argc; i++){
        if(argv[i][0]!='-')
            continue;
        switch(argv[i][1]){
            case 'n':
                n = atoi(argv[++i]);    // size of the dim 1 of the processor grid
                break;
            case 'm':
                m = atoi(argv[++i]);    // size of the dim 2 of the processor grid
                break; 
            case 'o':  
#ifdef _OPENMP
                o = atoi(argv[++i]);    // number of OpenMP threads
                omp_set_num_threads(o);
#endif
                break;
        }
    }

    parallel.initialize(n, m);

    const int sepLen=60;
    string SEP = "\n", solidSEP="\n"; 
    for(int i; i<sepLen; i++){
        SEP += "-";
        solidSEP += "_";
    }
    SEP += "\n\n"; 
    solidSEP += "\n\n";

    COUT << solidSEP;

    COUT << "Using " << n << " x " << m << " MPI processes; each spawning " << o << " OpenMP threads\n" << SEP;

    ///////////////////////////////////////////////////////////

    /* Set parameters that are not to be changed */

    const int dim       = 3;    // number of spatial dimensions
    const int halo      = 2;    // size of halo
    const int npts      = 500;  // number of lattice points in each direction
    const int nparts    = 400;  // number of particles in box
    const Real lat_res  = 0.1;  // lattice resolution 


    COUT << "Simulation parameters\n";
    
    COUT<< "> number of dimensions:     " << dim << "\n"
        << "> number of lattice points: " << npts << "\n"
        << "> size of halo:             " << halo << "\n"
        << "> number of particles:      " << nparts << "\n"
        << "> lattice resolution:       " << lat_res << "\n";

    COUT << solidSEP;

    ///////////////////////////////////////////////////////////

    MPI_timer TIMER(8);


    int loop_tag[]  = {0};      // for_each
    int proj_tag[]  = {1,2,3};  // scalar, vector, sym. tensor
    int fft_tag[]   = {4,5};    // forward, backward
    int part_tag[]  = {6,7};    // update vel., move particle    



#ifdef LOOPCORR
    loopCorrectionSimple(&TIMER, loop_tag, npts, halo);
#endif

#ifdef PROJECTION
    partMeshProjSimple(&TIMER, proj_tag, npts, halo, nparts, lat_res, dim);
#endif

#ifdef FOURIER
    fasterFourierTransformSimple(&TIMER, fft_tag, npts, halo, nparts, lat_res, dim);
#endif

#ifdef PART
    particleUpdateSimple(&TIMER, part_tag, npts, halo, nparts, lat_res, dim);
#endif




    COUT << SEP;


    return 0;
}
