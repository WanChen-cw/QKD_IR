#ifndef CHECKVALUEEXTRACTOR_HPP_
#define CHECKVALUEEXTRACTOR_HPP_

#include "CRC.hpp"
#include "FrozenSiteGenerator.hpp"

template <typename B = int>
class CheckValueExtractor {
private:
    
protected:
    std::shared_ptr<CRC<B>> crc;
public:
    CheckValueExtractor(const std::shared_ptr<CRC<B>>& crc);// 

    virtual ~CheckValueExtractor() = default;

    virtual void cvextract(const std::vector<B>& U, std::vector<B>& checkvalue) = 0;

    uint32_t crcextract(const std::vector<B>& U);

    std::vector<uint32_t> segcrcextract(const std::vector<B>& U, std::vector<int> &SegSite);
};

template<typename B>
CheckValueExtractor<B>
::CheckValueExtractor(const std::shared_ptr<CRC<B>>& crc)
    : crc(crc) {
}

template<typename B>
uint32_t CheckValueExtractor<B>
::crcextract(const std::vector<B>& U){
    return crc->calculate(U);
}

template<typename B>
std::vector<uint32_t> CheckValueExtractor<B>
::segcrcextract(const std::vector<B>& U, std::vector<int>& SegSite){
    return crc->segcrccalculate(U, SegSite);
}


template <typename B = int>
class CheckValueExtractorPolar :public CheckValueExtractor <B>{
private:

protected:
    std::shared_ptr<FrozenSiteGenerator> fbg;
    int N;
public:
    CheckValueExtractorPolar(const std::shared_ptr<CRC<B>>& crc,
        const std::shared_ptr<FrozenSiteGenerator>& fbg,
        int N)
        : CheckValueExtractor<B>(crc), fbg(fbg), N(N) {}

    virtual~CheckValueExtractorPolar() override {}

    virtual void cvextract(const std::vector<B>& U, std::vector<B>& checkvalue) override;
};

template<typename B>
void CheckValueExtractorPolar<B>
::cvextract(const std::vector<B>& U, std::vector<B>& checkvalue){
    for (int i = 0; i < N; ++i) {
        if (fbg->frozen_bits[i]) {
            checkvalue[i] = U[i];
        }
    }
}


template <typename B = int>
class CheckValueExtractorPcPolar :public CheckValueExtractorPolar <B> {
private:

protected:
    int p;
public:
    CheckValueExtractorPcPolar(const std::shared_ptr<CRC<B>>& crc,
        const std::shared_ptr<FrozenSiteGenerator>& fbg,
        int N,int p=9)
        : CheckValueExtractorPolar<B>(crc ,fbg, N),p(p) {}

    virtual ~CheckValueExtractorPcPolar() override {}

    virtual void cvextract(const std::vector<B>& U, std::vector<B>& checkvalue) override;
};

template<typename B>
void CheckValueExtractorPcPolar<B>
::cvextract(const std::vector<B>& U, std::vector<B>& checkvalue) {
    std::copy(U.begin(), U.begin()+p, checkvalue.begin());
    for (int i = p; i < this->N; ++i) {
        checkvalue[i] = checkvalue[i-p] ^ U[i];
    }
}


template <typename B = int>
class CheckValueExtractorPcPolarOpt :public CheckValueExtractorPcPolar <B> {
private:

protected:
    std::vector<B>   reg;
public:
    CheckValueExtractorPcPolarOpt(const std::shared_ptr<CRC<B>>& crc,
        const std::shared_ptr<FrozenSiteGenerator>& fbg,
        int N, int p = 9)
        : CheckValueExtractorPcPolar<B>(crc,fbg,N,p),reg(this->p,0) {}

    virtual ~CheckValueExtractorPcPolarOpt() override {}

    virtual void cvextract(const std::vector<B>& U, std::vector<B>& checkvalue) override;
};

template<typename B>
void CheckValueExtractorPcPolarOpt<B>
::cvextract(const std::vector<B>& U, std::vector<B>& checkvalue) {
    int idx = 0;
    std::fill(reg.begin(), reg.end(),0);
    for (int i = 0; i < this->N; ++i) {
        if (this->fbg->frozen_bits[i]) {
            checkvalue[i] = reg[idx] ^ U[i]; // XOR operation
        }
        else {
            reg[idx] = reg[idx] ^ U[i]; // XOR operation
            checkvalue[i] = reg[idx];
        }
        idx++;
        if (idx == this->p) {idx = 0; }
    }
}

#endif /*CHECKVALUEEXTRACTOR_HPP_ */


