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



void compactSaveField(string key, Field<Real>* field){
    int lat_size[] = {1, field->lattice().size(1), field->lattice().size(2)};
    Lattice dummy_lat(3, lat_size, field->lattice().halo());
    
    Field<Real> dummy_field(dummy_lat, 1);




    // #pragma omp parallel for collapse(2)
    for(int k=0; k<field->lattice().sizeLocal(2); k++)
        for(int j=0; j<field->lattice().sizeLocal(1); j++)
        {
            int ijk[] = {0,j,k};
            int jk[] = {j,k};

            Site x_dummy(dummy_lat, dummy_lat.indexTransform(jk));

            Site x(field->lattice(), field->lattice().indexTransform(ijk));
            Real s = 0.;
            
            for(int i=0; i< field->lattice().sizeLocal(0); i++){
                s += (*field)(x);
                x.indexAdvance(1);
            }
            dummy_field(x_dummy) = s;
        }
    

    
    int w = 1;
    string filename = getOutputFilename(key);
    dummy_field.saveHDF5(filename);//, 0, 1);

    COUT << "File '" << filename << "' was saved.\n";

    dummy_field.dealloc();


}


// void saveAll(Field<Real>* scalar_field, string key){


//     string filename = getOutputFilename(key);
//     scalar_field.saveHDF5(filename);

// }