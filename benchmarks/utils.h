#include <stdlib.h>
#include "LATfield2.hpp"
using namespace LATfield2;



class Diagnostics{
    private:

        ofstream epicrisis;
        string outfile = "epicrisis.txt";

        // information:
        string key;

        // key properties:
        int n, m;


        // performance metrics:
        double org_runtime;
        double new_runtime;
        double performance_gain;

        // error diagnostics:
        double max_error;
        double min_error;
        double avg_error;

    public:

        // constructors:

        Diagnostics(string key, int n, int m);

        void write_epicrisis();

        void error_diagnosis(double max, double min, double avg);

        void performance_metrics(double org_time, double new_time);

};


Diagnostics::Diagnostics(string key, int n, int m) : 
key(key), n(n), m(m){}


void Diagnostics::write_epicrisis(){
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






void field_comparison(Lattice lat, Diagnostics *d){

    parallel.barrier();

    string outfile = "fresh_output.h5";
    string orgfile = "org_output.h5";

    /* Read files and create fields */

    Field<Real> B_fresh;
    B_fresh.initialize(lat); 
    B_fresh.alloc();
    B_fresh.loadHDF5(outfile);

    Field<Real> B_org;
    B_org.initialize(lat); 
    B_org.alloc();
    B_org.loadHDF5(orgfile);

    /* Compare fields 
    (going with absolute error instead of relative error to avoid division by zero) */

    Site x(lat);

    double abs_err;         // absolute error
    double max_err = 0.;    // maximal error
    double min_err = 20.;   // minimal error
    double avg_err = 0.;    // average error

    for(x.first(); x.test(); x.next()){
        abs_err = fabs( B_fresh(x) - B_org(x) );
        avg_err += abs_err;
        if(min_err > abs_err) min_err = abs_err;
        if(max_err < abs_err) max_err = abs_err;   
    }

    parallel.max(max_err);
    parallel.min(min_err);
    parallel.sum(avg_err);

    avg_err /= (double)lat.sites();

    parallel.barrier();

    d->error_diagnosis(max_err, min_err, avg_err);
}