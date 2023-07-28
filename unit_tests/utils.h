#include <random>
#include "LATfield2.hpp"
using namespace LATfield2;

#include "../../LATfield2/timer.hpp"

string getOutputFilename(string key){
    string file = "output/";
#ifdef _OPENMP
    file += "new";
#else  
    file += "org";
#endif
    file += "Result_";
    file += key;
    file += ".h5";

    return file;

}




template <typename part, typename part_info, typename part_dataType>
void randomParticleEnsemble(Particles<part,part_info,part_dataType>* parts, const int nparts, Real *boxSize, unsigned int seed){

    srand(seed);
    Real norm = 1./(Real)RAND_MAX;
    for(long idx=0; idx<nparts; idx++){
        part p;
        p.ID = idx;
        for(int i=0; i<3; i++){
            p.pos[i] = (Real)static_cast<double>(rand()) * norm * (Real)boxSize[i];
            p.vel[i] = (Real)static_cast<double>(rand()) * pow(norm, 1.5) * (Real)boxSize[i];
        }

        parts->addParticle_global(p);
    }


}