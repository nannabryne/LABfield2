

void fasterFourierTransformSimple(MPI_timer *timer, int* run_no, const int npts=1000, const int halo=2, const int nparts=100, const Real lat_res=0.1, const int dim=3){


    Lattice lat_field_x(dim, npts, halo);
    Lattice lat_field_k;
    lat_field_k.initializeRealFFT(lat_field_x, 0);


    Field<Real> scalar_field_x(lat_field_x, 1);
    Field<Imag> scalar_field_k(lat_field_k, 1);


    double sigma_sqrd   = lat_field_x.size(0)*lat_field_x.size(1)/10; 
    double mu[3]        = {lat_field_x.size(0)*0.5, lat_field_x.size(1)*0.5, lat_field_x.size(2)*0.5};
    double A            = pow(2.*M_PI*sigma_sqrd, -dim*0.5); 

    auto gaussian = [&](Site x){
        double xp1, xp2, xp3;   // x' = x - μ = (x'_1, x'_2, x'_3) 
        xp1 = (x.coord(0) - mu[0]); xp1 *= xp1;
        xp2 = (x.coord(1) - mu[1]); xp2 *= xp2;
        xp3 = (x.coord(2) - mu[2]); xp3 *= xp3;

        scalar_field_x(x) = A*exp(-(xp1+xp2+xp3) / sigma_sqrd*.5);
    };

    lat_field_x.for_each(gaussian);

    /* Fast fourier transform */
    PlanFFT<Imag> plan(&scalar_field_x, &scalar_field_k);


    timer->start(run_no[0]);
    plan.execute(FFT_FORWARD);
    timer->stop(run_no[0]);



    timer->start(run_no[1]);
    plan.execute(FFT_BACKWARD);
    timer->stop(run_no[1]);

    scalar_field_x.updateHalo();
    




    /* Save result in compact manner */

    // Field<Real> dummy_field(lat_field_x, 1);


    int xcoord = npts/4;
    int w = 3;
    string filename = getOutputFilename("FasterFourierTransform");
    scalar_field_x.saveSliceHDF5(filename, xcoord, w);

    COUT << "File '" << filename << "' was saved.\n";


}