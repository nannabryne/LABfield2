#include <stdlib.h>
#include "LATfield2.hpp"
using namespace LATfield2;



class Diagnostics{
    private:

        ofstream epicrisis;
        string outfile = "epicrisis.txt";

        // information:
        string key;     // key descibing unit test

        // key properties:
        int n, m;           // parallel grid size
        int num_prcs;       // total number of processes
        int dim = 3;
        int num_cells;      // 
        int num_parts = 0;  // 

        // reference pararameters:
        int REF_num_prcs;
        double REF_runtime;

        // performance metrics:
        double runtime;
        double org_runtime;         // 'original' runtime
        double new_runtime;         // current runtime
        double performance_gain;    // gain in performance in % (loss if negative)
        double efficiency;

        // error diagnostics:
        double max_error;   // maximal error (absolute)
        double min_error;   // minimal error (absolute)
        double avg_error;   // average error (absolute)

        string text_epicrisis();    // ???


    public:

        // constructors:

        Diagnostics(string key, int n, int m);

        void write_epicrisis();

        void print_epicrisis();

        void error_diagnosis(double max, double min, double avg);

        void performance_metrics(double org_time, double new_time);

        void compute_performance_metrics(double comp_duration);

        void provide_reference_parameters(int num_prcs, double comp_duration);

        void provide_computation_parameters(int num_cells_per_dim, int num_parts_per_dim=0, int num_dimensions=3);

};


Diagnostics::Diagnostics(string key, int n, int m) : 
key(key), n(n), m(m){
    num_prcs = n*m;
}

void Diagnostics::provide_reference_parameters(int num_prcs, double comp_duration){
    REF_num_prcs = num_prcs;
    REF_runtime = comp_duration;
}

void Diagnostics::provide_computation_parameters(int num_cells_per_dim, int num_parts_per_dim, int num_dimensions){
    num_cells = num_cells_per_dim;
    num_parts = num_parts_per_dim;
    dim = num_dimensions;
}


void Diagnostics::write_epicrisis(){
    if(parallel.isRoot()){
    epicrisis.open(outfile.c_str());

    epicrisis << "Basic information" << endl
                << "benchmark key: " << key << endl
                << "operations: " << "..." << endl 
                << endl;

    epicrisis << "Key parameters" << endl
                << "parallel grid: " << "(" << n << " x " << m << ")" << endl
                << "lattice size: " << num_cells << "^" << dim << endl
                << "number of particles: " << num_parts << "^" << dim << endl
                << endl;

    epicrisis << "Performance metrics" << endl
                << "original runtime: " << REF_runtime << " ms" << endl
                << "latest runtime: " << runtime << " ms" << endl
                << "performance enhancement: " << performance_gain << " %" << endl
                << "efficiency: " << efficiency << endl
                << endl;
    
    epicrisis << "Error diagnostics" << endl
                << "maximal error: " << max_error << endl
                << "minimal error: " << min_error << endl
                << "average error: " << avg_error << endl
                << endl;
        
    epicrisis.close();
    }
    
}


void Diagnostics::print_epicrisis(){ // temporary solution

    COUT << "*******************************" << endl << endl;

    COUT << "Basic information" << endl
                << "benchmark key: " << key << endl
                << "operations: " << "..." << endl 
                << endl;

    COUT << "Key parameters" << endl
                << "parallel grid: " << "(" << n << " x " << m << ")" << endl
                << "lattice size: " << num_cells << "^" << dim << endl
                << "number of particles: " << num_parts << "^" << dim << endl
                << endl;

    COUT << "Performance metrics" << endl
                << "original runtime: " << REF_runtime << " ms" << endl
                << "latest runtime: " << runtime << " ms" << endl
                << "performance enhancement: " << performance_gain << " %" << endl
                << "efficiency: " << efficiency << endl
                << endl;
    
    COUT << "Error diagnostics" << endl
                << "maximal error: " << max_error << endl
                << "minimal error: " << min_error << endl
                << "average error: " << avg_error << endl
                << endl;

    COUT << endl << "*******************************" << endl;
    
}

void Diagnostics::error_diagnosis(double max, double min, double avg){
    max_error = max;
    min_error = min;
    avg_error = avg;
}


void Diagnostics::compute_performance_metrics(double comp_duration){

    runtime = comp_duration;
    performance_gain = (REF_runtime-runtime)/REF_runtime*100.0;
    if(fabs(performance_gain) < 8.){
        performance_gain = 0.;
    }
    efficiency = (REF_num_prcs*REF_runtime)/(num_prcs*runtime);
  
}






void field_comparison(int dim, int lat_size, int halo_size, Diagnostics *d){

    // NB! Only works with scalar fields...

    string outfile = "fresh_output.h5";
    string orgfile = "org_output.h5";

    Lattice lat(dim, lat_size, halo_size);

    /* Read files and create fields */


    Field<Real> field_fresh(lat);
    field_fresh.loadHDF5(outfile);

    Field<Real> field_org(lat);
    field_org.loadHDF5(orgfile);

    /* Compare fields 
    (going with absolute error instead of relative error to avoid division by zero) */

    Site x(lat);

    double err;
    double abs_err;         // absolute error
    double rel_err;         // relative error
    double max_err = 0.;    // maximal error
    double min_err = 1.0e4; // minimal error
    double avg_err = 0.;    // average error

    for(x.first(); x.test(); x.next()){
        abs_err = fabs( field_fresh(x) - field_org(x) );
        rel_err = abs_err/fabs(field_org(x));
        err = abs_err;
        avg_err += err;
        if(min_err > err) min_err = err;
        if(max_err < err) max_err = err;   
    }

    parallel.max(max_err);
    parallel.min(min_err);
    parallel.sum(avg_err);

    avg_err /= (double)lat.sites();

    parallel.barrier();

    d->error_diagnosis(max_err, min_err, avg_err);
}