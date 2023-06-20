/**
 * @file fft_execution.hpp
 * @author Nanna Bryne
 * 
 * @brief Unit test: run code from '../boilerplate_examples/fft_execution.hpp' and get diagnostic data
 * 
*/





#include "../boilerplate_examples/fft_execution.hpp"


class Diagnostics{
    private:

        const string output_path = "../results/"; // for now
        string filename_ascii;
        ofstream epicrisis;

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

        Diagnostics(string filename_ascii, int n, int m);

        void write_epicrisis();

        void error_diagnosis(double max, double min, double avg);

        void performance_metrics(double org_time, double new_time);


};


Diagnostics::Diagnostics(string filename_ascii, int n, int m) : 
filename_ascii(filename_ascii), n(n), m(m){
    key = filename_ascii;
}

void Diagnostics::write_epicrisis(){
    const string filename = filename_ascii + ".txt";
    epicrisis.open(output_path + filename.c_str());

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



int main(int argc, char **argv){

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
  
    parallel.initialize(n, m);
    COUT << " " << endl;
    

    double runtime;
    run_example(&runtime);

    Diagnostics diag("test", n, m);
    

    diag.error_diagnosis(1., 0., .6);
    diag.performance_metrics(5.9, runtime);
    
    diag.write_epicrisis();
    
    
    // End timing!
    // double stop = MPI_Wtime();

   

    /* Compare results */




    /* Run diagnostics */

    // double duration = (stop-start)*1000.0;

    COUT << "runtime: "<< runtime << " ms"<< endl;
    // diagnostics.new_runtime = duration;



    COUT << " " << endl;




    return 0;
}



