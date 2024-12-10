#include <iostream>
#include <memory>
#include <vector>
#include <cmath>
#include "SiftKeyGenerator.hpp"
#include "Encoder.hpp"
#include "CRC.hpp"
#include "FrozenSiteGenerator.hpp"
#include "CheckValueExtractor.hpp"
#include "Decoder.hpp"
#include "Display.hpp"

#include <thread>
#include <vector>

class IR {
private:
    
public:
    IR(int N, int L, int seed, float qber, const std::string& fbgpath, int Nf, uint32_t crcPolynomial, std::string& type,Reporter& reporter);

    void run();
    virtual void informationrecongnize();

protected:
    //parament
    int N;
    int L;
    int seed;
    float qber;
    std::string fbgpath;
    int Nf;
    uint32_t crcPolynomial;
    const std::string type;

    std::unique_ptr<SiftKeyGeneratorContext<int>> siftkey;
    std::unique_ptr<Encoder<int>> encoder;
    std::shared_ptr<CRC<int>> crc;
    std::shared_ptr<FrozenSiteGenerator> fbg;
    std::unique_ptr<CheckValueExtractor<int>> extractor;
    std::unique_ptr<Decoder<int>> decoder;

    Reporter& reporter;

    //buffer
    std::vector<int> ksa;
    std::vector<int> ksb;
    std::vector<int> U;
    uint32_t crc32;
    std::vector<int> checkvalue;
    std::vector<int> _U;
    size_t testnum;


    virtual void init();
};

IR::IR(int N, int L, int seed, float qber, const std::string& fbgpath, int Nf, uint32_t crcPolynomial, std::string& type,Reporter& reporter)
    : N(N), L(L), seed(seed), qber(qber), fbgpath(fbgpath), Nf(Nf), crcPolynomial(crcPolynomial), reporter(reporter), type(type),
    ksa(N), ksb(N), U(N), crc32(0), checkvalue(N), _U(N), testnum(0){
    siftkey = std::make_unique<SiftKeyGeneratorContext<int>>();
    siftkey->setStrategy(std::make_unique<RandomSiftKeyGenerator<int>>(N, qber, seed));
    encoder = std::make_unique<EncoderPolar<int>>();
    crc = std::make_shared<CRC<int>>(crcPolynomial);
    fbg = std::make_shared<FrozenSiteGenerator>(N, fbgpath);

    
}

void IR::run() {
    this->init();
    while(reporter.iscontinue()){
        testnum++;
        siftkey->generate(ksa, ksb);
        encoder->encode(ksa, U);
        informationrecongnize();
    }
}

void IR
::informationrecongnize(){
    fbg->generate(Nf);
    crc32 = extractor->crcextract(U);
    extractor->cvextract(U, checkvalue);

    if (!decoder->decode(ksb, checkvalue, crc32, _U)) {
        reporter.FEAdd(testnum);
        testnum = 0;
    }
}

void IR
::init(){
    if (type == "PCPOLAR") {
        extractor = std::make_unique<CheckValueExtractorPcPolarOpt<int>>(crc, fbg, N);
        decoder = std::make_unique<DecoderPcPolar<int>>(N, L, qber, crc, fbg);
    }
    else if(type == "POLAR"){
        extractor = std::make_unique<CheckValueExtractorPolar<int>>(crc, fbg, N);
        decoder = std::make_unique<DecoderPolar<int>>(N, L, qber, crc, fbg);
    }
}


class AIR :public IR {
private:

public:
    AIR(int N, int L, int seed, float qber, const std::string& fbgpath, int Nf, uint32_t crcPolynomial, std::string& type, Reporter& reporter,int t);

    ~AIR();

    virtual void informationrecongnize()override;
protected:
    int t;
};

AIR::AIR(int N, int L, int seed, float qber, const std::string& fbgpath, int Nf, uint32_t crcPolynomial, std::string& type, Reporter& reporter,int t)
    : IR(N, L, seed, qber, fbgpath, Nf, crcPolynomial, type, reporter),t(t){

}

AIR
::~AIR(){
}

void AIR
::informationrecongnize(){
    fbg->generate(Nf);
    crc32 = extractor->crcextract(U);
    extractor->cvextract(U, checkvalue);

    int addrounds = 0;
    while (!decoder->decode(ksb, checkvalue, crc32, _U)) { 
        fbg->addfb(t);
        addrounds++;
        extractor->cvextract(U, checkvalue);
    }
    reporter.FNAdd(addrounds);
}


class SAIR :public AIR {
private:

public:
    SAIR(int N, int L, int seed, float qber, const std::string& fbgpath, int Nf, uint32_t crcPolynomial, std::string& type, Reporter& reporter, int t, int Ns, int Nf_max);

    ~SAIR();

    void informationrecongnize()override;
protected:
    int Ns;
    int Nf_max;
    std::vector<uint32_t> crc_seg;


    void init()override;
};

SAIR::SAIR(int N, int L, int seed, float qber, const std::string& fbgpath, int Nf, uint32_t crcPolynomial, std::string& type, Reporter& reporter, int t, int Ns,int Nf_max)
    : AIR(N, L, seed, qber, fbgpath, Nf, crcPolynomial, type,reporter, t),Ns(Ns), Nf_max(Nf_max), crc_seg(Ns,0){
    
    
}

SAIR
::~SAIR() {
}

void SAIR
::informationrecongnize() {
    fbg->generate(Nf);
    crc32 = extractor->crcextract(U);
    fbg->crc_seg = extractor->segcrcextract(U, fbg->SegSite);
    extractor->cvextract(U, checkvalue);

    int addrounds = 0;
    int Nf_t = 0;
    while (!decoder->decode(ksb, checkvalue, crc32, _U)) {///daigai 
        Nf_t+=fbg->addfb_s(t);//daigai
        addrounds++;
        extractor->cvextract(U, checkvalue);
    }
    reporter.FNsAdd(addrounds, Nf_t);//daigai
}

void SAIR
::init(){
    if (type == "PCPOLAR") {
        extractor = std::make_unique<CheckValueExtractorPcPolarOpt<int>>(crc, fbg, N);
        decoder = std::make_unique<DecoderPcPolar<int>>(N, L, qber, crc, fbg);
    }
    else if (type == "POLAR") {
        extractor = std::make_unique<CheckValueExtractorPolar<int>>(crc, fbg, N);
        decoder = std::make_unique<DecoderPolar<int>>(N, L, qber, crc, fbg);
    }
    else if (type == "SPCPOLAR") {
        extractor = std::make_unique<CheckValueExtractorPcPolarOpt<int>>(crc, fbg, N);
        decoder = std::make_unique<DecoderSegPcPolar<int>>(N, L, qber, crc, fbg);
        fbg->generate(this->Nf);
        this->fbg->SegSiteGenerate(this->Ns, this->Nf, this->Nf_max);
    }
}
