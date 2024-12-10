#include <omp.h> 
#include <iostream>

#include "IR.hpp"


int main() {
    const int numThreads = 20; 

    int N=65536;
    int L=16;
    int seed=12;
    float QBER=0.05;
    std::string fbgpath= "conf/N_65536_005.pc";
    int Nf=20350;
    uint32_t crcPolynomial= 0x04C11DB7;


    int t = 100;
    int Ns = 4;
    int Nf_max = 22500;


    Reporter reporter(N, L, seed, QBER, fbgpath, Nf, crcPolynomial, t);
    reporter.printparament();

#pragma omp parallel num_threads(numThreads)
    {

#pragma omp single
       {
           std::cout << "Using " << omp_get_num_threads() << " threads." << std::endl;
           reporter.set_num_threads(omp_get_num_threads());
       }
        //std::string type = "PCPOLAR";//POLAR PCPOLAR SPCPOLAR
        //IR test(N, L, seed, QBER, fbgpath, Nf, crcPolynomial,type, reporter);

        std::string type = "PCPOLAR";//POLAR PCPOLAR SPCPOLAR
        AIR test(N, L, seed, QBER, fbgpath, Nf, crcPolynomial, type, reporter,t);

        //std::string type = "SPCPOLAR";//POLAR PCPOLAR SPCPOLAR
        //SAIR test(N, L, seed, QBER, fbgpath, Nf, crcPolynomial, type,reporter, t, Ns, Nf_max);
        test.run();
    }


    return 0;
}
