

void partMeshProjSimple(MPI_timer *timer, int* run_no, const int npts=1000, const int halo=2, const int nparts=100, const Real lat_res=0.1, const int dim=3){


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

    unsigned int seed = 30;
    randomParticleEnsemble(&parts, nparts, boxSize, seed);


    /* Create fields */
    Field<Real> scalar_field(lat_field, 1);
    Field<Real> vector_field(lat_field, 3);
    Field<Real> tensor_field(lat_field, 3, 3, LATfield2::symmetric);    


    /* Projections */

    projection_init(&scalar_field);
    projection_init(&vector_field);
    projection_init(&tensor_field);


    timer->start(run_no[0]);
    scalarProjectionCIC_project(&parts, &scalar_field);
    timer->stop(run_no[0]);
    scalarProjectionCIC_comm(&scalar_field);

    timer->start(run_no[1]);
    vectorProjectionCICNGP_project(&parts, &vector_field);
    timer->stop(run_no[1]);
    vectorProjectionCICNGP_comm(&vector_field);

    timer->start(run_no[2]);
    symtensorProjectionCICNGP_project(&parts, &tensor_field);
    timer->stop(run_no[2]);
    symtensorProjectionCICNGP_comm(&tensor_field);





    
    
    /* Save result in compact manner */


    Field<Real> dummy_field(lat_field, 1);

    auto combine = [&](Site& x){
        Real s = scalar_field(x);
        Real v = vector_field(x,0) + vector_field(x,1) + vector_field(x,2);
        Real t = tensor_field(x,0,0) + tensor_field(x,1,1);
        dummy_field(x) = s + v + t;
    };
    lat_field.for_each(combine);

    int xcoord = npts/4;
    int w = 3;
    string filename = getOutputFilename("PartMeshProjection");
    scalar_field.saveSliceHDF5(filename, xcoord, w);

    COUT << "File '" << filename << "' was saved.\n";






}