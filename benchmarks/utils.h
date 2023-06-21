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
        int n, m;       // parallel grid size

        // performance metrics:
        double org_runtime;         // 'original' runtime
        double new_runtime;         // current runtime
        double performance_gain;    // gain in performance in % (loss if negative)

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

};


Diagnostics::Diagnostics(string key, int n, int m) : 
key(key), n(n), m(m){}



void Diagnostics::write_epicrisis(){
    if(parallel.isRoot()){
    epicrisis.open(outfile.c_str());

    epicrisis << "Basic information" << endl
                << "benchmark key: " << key << endl
                << "operations: " << "..." << endl 
                << endl;

    epicrisis << "Key parameters" << endl
                << "parallel grid: " << "(" << n << ", " << m << ")" << endl
                << "lattice size: " << "..." << endl
                << endl;

    epicrisis << "Performance metrics" << endl
                << "original runtime: " << org_runtime << " ms" << endl
                << "latest runtime:   " << new_runtime << " ms" << endl
                << "performance enhancement: " << performance_gain << " %" << endl
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
                << "parallel grid: " << "(" << n << ", " << m << ")" << endl
                << "lattice size: " << "..." << endl
                << endl;

    COUT << "Performance metrics" << endl
                << "original runtime: " << org_runtime << " ms" << endl
                << "latest runtime:   " << new_runtime << " ms" << endl
                << "performance enhancement: " << performance_gain << " %" << endl
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

void Diagnostics::performance_metrics(double org_time, double new_time){
    org_runtime = org_time;
    new_runtime = new_time;
    performance_gain = (org_runtime-new_runtime)/org_runtime*100.0;
    if(fabs(performance_gain) < 5.){
        performance_gain = 0.;
    }
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