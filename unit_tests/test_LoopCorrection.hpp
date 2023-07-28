


void loopCorrectionSimple(MPI_timer *timer, int* run_no, const int npts=1000, const int halo=2, const int nparts=100, const Real lat_res=0.1, const int dim=3){


    Lattice lat_field(dim, npts, halo);


    Field<Real> scalar_field(lat_field, 1);

    auto op = [&](Site &x){
        scalar_field(x) = x.index();
    };

    timer->start(run_no[0]);
    lat_field.for_each(op);
    timer->stop(run_no[0]);


    int xcoord = npts/4;
    int w = 3;
    string filename = getOutputFilename("LoopCorrection");
    scalar_field.saveSliceHDF5(filename, xcoord, w);

    COUT << "File '" << filename << "' was saved.\n";


}