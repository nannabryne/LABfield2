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
    Diagnostics d("particle_mesh_projection", n, m);

    const int dim         = 3;      // #{dimensions}      
    const int lat_size    = 60;     // size of lattice in one dim.  FIXME
    const int npts        = 60;     // size of lattice in one dim. 
    const int halo_size   = 1;      // size of halo in lattice
    
    
    int num_parts   = 60;       // #{particles}   
    Real lat_res    = 0.1;      // lattice resolution
    
    
    
    parallel.initialize(n, m);

    COUT << " " << endl;


    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


    /* Initialise lattices and fields */

    Lattice lat(dim, lat_size, halo_size);
    Lattice lat_part(dim, lat_size, 0);

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
    int ratio = num_parts/lat_size;

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

    /* Projection initialisation */

    projection_init(&phi);
    projection_init(&B);
    projection_init(&T);


    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

    /* Actual projection*/

    // ! START TIMING HERE
    timer_start = MPI_Wtime();

    scalarProjectionCIC_project(&parts, &phi);
    scalarProjectionCIC_comm(&phi);

    vectorProjectionCIC_project(&parts, &B);
    vectorProjectionCIC_comm(&B);

    symtensorProjectionCICNGP_project(&parts, &T);
    symtensorProjectionCICNGP_comm(&T);


    // ! END TIMING HERE
    runtime = ((MPI_Wtime() - timer_start)*1000.0);


    phi.saveHDF5(outfile);





    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

    /* Ensure validity of output output */

    parallel.barrier();
    field_comparison(dim, lat_size, halo_size, &d);

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

    /* Diagnostic analysis*/

    d.performance_metrics(17, runtime);

    d.write_epicrisis();
    d.print_epicrisis();

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


    COUT << " " << endl;

    return 0;
}

