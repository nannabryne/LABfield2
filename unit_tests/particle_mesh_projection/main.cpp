/**
 * @file particle_mesh_projection/main.cpp
 * @author Nanna Bryne
 * 
 * @brief Unit test: ...
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
    const int dim         = 3;      // #{dimensions}      
    const int npts        = 512;    // size of lattice in one dim. 
    const int halo        = 1;      // size of halo in lattice
    const int num_parts   = 512;    // #{particles}   

    Real lat_res    = 0.1;       // lattice resolution
    
    Diagnostics d("particle_mesh_projection", n, m);
    d.provide_reference_parameters(64, 4000);
    d.provide_computation_parameters(npts, num_parts);
    
    
    parallel.initialize(n, m);

    COUT << " " << endl;


    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


    /* Initialise lattices and fields */

    Lattice lat(dim, npts, halo);
    Lattice lat_part(dim, npts, 0);

    Field<Real> phi(lat, 1);                            // scalar field φ
    Field<Real> B(lat, dim);                            // vector field B
    Field<Real> T(lat, dim, dim, LATfield2::symmetric); // tensor field T

    Real box_size[3];
    for(int i=0; i<3; i++)box_size[i] = lat_res * lat_part.size(i);



    Site x(lat);

    for(x.first(); x.test(); x.next()){
        phi(x) = 0.;
        for(int i=0; i<3; i++)B(x,i)=0;
    }



    /* Create particles 
        (basically replica of 'particles.cpp') */

    part_simple_info particles_global_info;
    part_simple_dataType particles_dataType;

    particles_global_info.mass = 0.2;
    particles_global_info.relativistic = false;
    set_parts_typename(&particles_global_info, "part_simple");

    Particles<part_simple,part_simple_info,part_simple_dataType> parts;
    parts.initialize(particles_global_info, particles_dataType, &lat_part, box_size);

    part_simple part;

    long index = 0;
    int ratio = num_parts/npts;

    Site x_part(lat_part);
    for(x_part.first(); x_part.test(); x_part.next()){

        for(int i=0; i<ratio; i++)
            for(int j=0; j<ratio; j++)
                for(int k=0; k<ratio; k++){

                    part.ID = index;

                    part.pos[0] = (Real)x_part.coord(0) * (Real)box_size[0] / (Real)npts;
                    part.pos[1] = (Real)x_part.coord(1) * (Real)box_size[1] / (Real)npts;
                    part.pos[2] = (Real)x_part.coord(2) * (Real)box_size[2] / (Real)npts;

                    part.vel[0] = 1.;
                    part.vel[1] = 1.;
                    part.vel[2] = 1.;

                    parts.addParticle_global(part);
                    index++;

                }
    }

    /* Projection initialisation
        (set the field to zero everyhvere) */

    projection_init(&phi);
    projection_init(&B);
    projection_init(&T);


    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

    /* Actual projection */

    double timer_pause, timer_play;

    // ! START TIMING HERE
    timer_start = MPI_Wtime();

    COUT << endl << "scalar projection ..." << endl;
    timer_play = MPI_Wtime();
    scalarProjectionCIC_project(&parts, &phi);
    timer_pause = MPI_Wtime();
    COUT << "finished in " << (timer_pause-timer_play)*1000.0 << " ms" << endl;
    
    COUT << endl << "scalar communication ..." << endl;
    timer_play = MPI_Wtime();
    scalarProjectionCIC_comm(&phi);
    timer_pause = MPI_Wtime();
    COUT << "finished in " << (timer_pause-timer_play)*1000.0 << " ms" << endl;




    COUT << endl << "vector projection ..." << endl;
    timer_play = MPI_Wtime();
    vectorProjectionCIC_project(&parts, &B);
    timer_pause = MPI_Wtime();
    COUT << "finished in " << (timer_pause-timer_play)*1000.0 << " ms" << endl;
    
    COUT << endl << "vector communication ..." << endl;
    timer_play = MPI_Wtime();
    vectorProjectionCIC_comm(&B);
    timer_pause = MPI_Wtime();
    COUT << "finished in " << (timer_pause-timer_play)*1000.0 << " ms" << endl;




    COUT << endl << "symmetric tensor projection ..." << endl;
    timer_play = MPI_Wtime();
    symtensorProjectionCICNGP_project(&parts, &T);
    timer_pause = MPI_Wtime();
    COUT << "finished in " << (timer_pause-timer_play)*1000.0 << " ms" << endl;
    
    COUT << endl << "symmetric tensor communication ..." << endl;
    timer_play = MPI_Wtime();
    symtensorProjectionCICNGP_comm(&T);
    timer_pause = MPI_Wtime();
    COUT << "finished in " << (timer_pause-timer_play)*1000.0 << " ms" << endl;



    // ! END TIMING HERE
    runtime = ((MPI_Wtime() - timer_start)*1000.0);


    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

    bool IS_REFERENCE = false;
    if(IS_REFERENCE){
        phi.saveHDF5(orgfile);
        COUT << endl << "reference runtime: " << runtime << " ms" << endl;
        COUT << "used " << n*m << " processes" << endl;
    }
    else{
    
    phi.saveHDF5(outfile);

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

    /* Ensure validity of output output */

    parallel.barrier();
    field_comparison(dim, npts, halo, &d);

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

    /* Diagnostic analysis*/


    d.compute_performance_metrics(runtime);

    d.write_epicrisis();
    d.print_epicrisis();

    }
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
    
   

    COUT << " " << endl;

    return 0;
}

