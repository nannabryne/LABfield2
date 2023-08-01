
Real update_velocity_func(  
    double dtau,
    double lat_resolution,
    part_simple * part,
    double * ref_dist,
    part_simple_info partInfo,
    Field<Real> ** fields,
    Site * sites,
    int nfields,
    double * params,
    double * outputs,
    int noutputs);

void move_particles_func(   
    double dtau,
    double lat_resolution,
    part_simple * part,
    double * ref_dist,
    part_simple_info partInfo,
    Field<Real> ** fields,
    Site * sites,
    int nfields,
    double * params,
    double * outputs,
    int noutputs);


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
    

    parts.updateVel(&update_velocity_func, 1.0);


    parts.moveParticles(&move_particles_func, 1.0);


    /* Save result in simple & compact way */

    Field<Real> dummy_field(lat_field);
    
    projection_init(&dummy_field);

    scalarProjectionCIC_project(&parts, &dummy_field);
    scalarProjectionCIC_comm(&dummy_field);

    int xcoord = npts/2;
    int w = 2*npts/3;
    string filename = getOutputFilename("ParticleUpdate");
    dummy_field.saveSliceHDF5(filename, xcoord, w);

    COUT << "File '" << filename << "' was saved.\n";

}   





Real update_velocity_func(  
    double dtau,
    double lat_resolution,
    part_simple * part,
    double * ref_dist,
    part_simple_info partInfo,
    Field<Real> ** fields,
    Site * sites,
    int nfields,
    double * params,
    double * outputs,
    int noutputs){
                    

    Real v2;

    for(int i=0; i<3; i++){
        Real f=0;
        for(int n=0; n<nfields; n++)f += (*fields[n])(sites[n], 0);
        (*part).vel[i] = (*part).vel[i] * (1.+f);

        v2 =(*part).vel[i] * (*part).vel[i];
    }

    return v2;
}
                          

void move_particles_func(   
    double dtau,
    double lat_resolution,
    part_simple * part,
    double * ref_dist,
    part_simple_info partInfo,
    Field<Real> ** fields,
    Site * sites,
    int nfields,
    double * params,
    double * outputs,
    int noutputs){

    
    for(int l=0; l<3; l++){
        Real f=0;
        for(int n=0; n<nfields; n++)f += (*fields[n])(sites[n], 0);
        (*part).pos[l] += dtau*(*part).vel[l]* (1. + 1e-4*f);
    }



}