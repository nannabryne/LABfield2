
void particleUpdateSimple(MPI_timer *timer, int* run_no, const int npts=1000, const int halo=2, const int nparts=100, const Real lat_res=0.1, const int dim=3){
    
    Real boxSize[dim];
    for(int i=0; i<dim; i++)boxSize[i] = lat_res*npts;

    Lattice lat_field(dim, npts, halo);
    Lattice lat_parts(dim, npts, 0);

    part_simple_info particles_global_info;
    part_simple_dataType particles_dataType;

    particles_global_info.mass = 0.2;
    particles_global_info.relativistic = false;
    set_parts_typename(&particles_global_info, "part_simple");

    Particles<part_simple,part_simple_info,part_simple_dataType> parts;
    parts.initialize(particles_global_info, particles_dataType, &lat_parts, boxSize);

    unsigned int seed = 140;
    randomParticleEnsemble(&parts, nparts, boxSize, seed);


    

}   