#ifndef DECODER_HPP_
#define DECODER_HPP_

#include <map>
#include <iostream>
#include <vector>
#include <stack>
#include <cmath>
#include <algorithm>
#include <numeric>
#include "CRC.hpp"
#include "FrozenSiteGenerator.hpp"

template <typename B = int>
class Decoder {
public:
    Decoder() {}

    virtual ~Decoder() = default;

    virtual bool decode(const std::vector<B>& ksb, const std::vector<B>& checkvalue, uint32_t crc32, std::vector<B>& _U) = 0;
};


template <typename B = int, typename R = float>
class DecoderPolar: public Decoder<B> {
private:
    
public:
    DecoderPolar(const int& N, const int& L, const float QBER, const std::shared_ptr<CRC<B>>& crc,const std::shared_ptr<FrozenSiteGenerator>& fbg);

    ~DecoderPolar();

    bool decode(const std::vector<B>& ksb, const std::vector<B>& checkvalue, uint32_t crc32, std::vector<B>& _U) override;
protected:
    //parament and tool
    int N;
    int m;
    int L;
    float QBER;
    float llr0;
    float llr1;
    std::shared_ptr<FrozenSiteGenerator> fbg;
    std::shared_ptr<CRC<B>> crc;
    std::vector<int> lambda_offset;
    std::vector<int> llr_layer_vec;
    std::vector<int> bit_layer_vec;

    //function
    void get_lambdaoffset();
    void get_bit_layer();
    void get_llr_layer();
    R f(R L1, R L2);
    void bittollr(const std::vector<B>& U, std::vector<R>& llr);
    virtual void decodeinit() {}
    virtual void unfrozenbitoperation(int phi, int phi_mod_2);
    virtual void frozenbitoperation(int phi ,int phi_mod_2,const std::vector<B>& checkvalue);
    //buffer
    std::vector<R>   llr;
    std::vector<std::vector<int>>       lazy_copy;
    std::vector<std::vector<R>>         P;
    std::vector<std::vector<B>>         C;
    std::vector<std::vector<B>>         u;
    std::vector<R>                      PM;
    std::vector<int>                    activepath;
    std::vector<R>                      PM_pair;
    std::vector<size_t>                 PM_inx;
    std::vector<int>                    compare;
};

template <typename B, typename R>
DecoderPolar<B, R>
::DecoderPolar(const int& N, const int& L, const float QBER, const std::shared_ptr<CRC<B>>& crc, const std::shared_ptr<FrozenSiteGenerator>& fbg)
    : N(N),L(L), QBER(QBER), fbg(fbg),crc(crc),m(log2(N)),bit_layer_vec(N, 0), llr_layer_vec(N, 0), llr(N, 0.0) , lazy_copy(m, std::vector<int>(L, 0)),
    P(N - 1, std::vector<R>(L, 0.0)) , C(N - 1, std::vector<B>(2 * L, 0)), u(L, std::vector<B>(N, 0)) , PM(L, 0.0), activepath(L, 0),
    PM_pair(2 * L), PM_inx(2 * L), compare(2 * L) {
    llr0 = log2((1 - QBER) / QBER);
    llr1 = log2(QBER / (1 - QBER));
    this->get_lambdaoffset();
    this->get_bit_layer();
    this->get_llr_layer();
}

template <typename B, typename R>
DecoderPolar<B, R>
::~DecoderPolar(){
}

template<typename B, typename R>
bool DecoderPolar<B, R>
::decode(const std::vector<B>& ksb, const std::vector<B>& checkvalue, uint32_t crc32, std::vector<B>& _U){
    bool IsSucessful = false;
    bittollr(ksb, llr);
    std::fill(PM.begin(), PM.end(), 0);
    std::fill(activepath.begin(), activepath.end(), 0);
    activepath[0] = 1;
    for (int i = 0; i < m; ++i) {
        lazy_copy[i][0] = 0;
    }
    this->decodeinit();
    int layer;
    int phi_mod_2;
    int index_1;
    int index_2;

    for (int phi = 0; phi < N; ++phi) {
        layer = llr_layer_vec[phi];
        phi_mod_2 = phi % 2;
        for (int l_index = 0; l_index < L; ++l_index) {
            if (activepath[l_index] == 0)
                continue;
            if (phi == 0) {
                index_1 = lambda_offset[m - 1];
                for (int beta = 0; beta < index_1; ++beta) {
                    P[beta + index_1 - 1][l_index] = f(llr[beta], llr[beta + index_1]);
                }

                for (int i_layer = m - 2; i_layer >= 0; --i_layer) {
                    index_1 = lambda_offset[i_layer];
                    index_2 = lambda_offset[i_layer + 1];
                    for (int beta = 0; beta < index_1; ++beta) {
                        P[beta + index_1 - 1][l_index] = f(P[beta + index_2 - 1][l_index], P[beta + index_1 + index_2 - 1][l_index]);
                    }
                }
            }
            else if (phi == N / 2) {
                index_1 = lambda_offset[m - 1];
                for (int beta = 0; beta < index_1; ++beta) {
                    int x_tmp = C[beta + index_1 - 1][2 * l_index];
                    P[beta + index_1 - 1][l_index] = (1 - 2 * x_tmp) * llr[beta] + llr[beta + index_1];
                }
                for (int i_layer = m - 2; i_layer >= 0; --i_layer) {
                    index_1 = lambda_offset[i_layer];
                    index_2 = lambda_offset[i_layer + 1];
                    for (int beta = 0; beta < index_1; ++beta) {
                        P[beta + index_1 - 1][l_index] = f(P[beta + index_2 - 1][l_index], P[beta + index_1 + index_2 - 1][l_index]);
                    }
                }

            }
            else {//----------------------
                index_1 = lambda_offset[layer];
                index_2 = lambda_offset[layer + 1];
                for (int beta = 0; beta < index_1; ++beta) {
                    P[beta + index_1 - 1][l_index] = (1 - 2 * C[beta + index_1 - 1][2 * l_index]) *
                        P[beta + index_2 - 1][lazy_copy[layer + 1][l_index]] + P[beta + index_1 + index_2 - 1][lazy_copy[layer + 1][l_index]];
                }
                for (int i_layer = layer - 1; i_layer >= 0; --i_layer) {
                    index_1 = lambda_offset[i_layer];
                    index_2 = lambda_offset[i_layer + 1];
                    for (int beta = 0; beta < index_1; ++beta) {
                        P[beta + index_1 - 1][l_index] = f(P[beta + index_2 - 1][l_index], P[beta + index_1 + index_2 - 1][l_index]);
                    }
                }
            }
        }

        //if now we decode an unfrozen bit
        if (!fbg->frozen_bits[phi]) {
            unfrozenbitoperation(phi, phi_mod_2);
        }
        //frozen bit operation
        else {
            frozenbitoperation(phi, phi_mod_2,checkvalue);
        }

        //partial-sum return
        for (int l_index = 0; l_index < L; ++l_index) {
            if (activepath[l_index] == 0)
                continue;
            if ((phi_mod_2 == 1) && (phi != N - 1)) {
                layer = bit_layer_vec[phi];
                for (int i_layer = 0; i_layer < layer; ++i_layer) {
                    index_1 = lambda_offset[i_layer];
                    index_2 = lambda_offset[i_layer + 1];
                    for (int beta = index_1; beta < 2 * index_1; ++beta) {
                        C[beta + index_1 - 1][2 * l_index + 1] = (C[beta - 1][2 * lazy_copy[i_layer][l_index]] + C[beta - 1][2 * l_index + 1]) % 2; // ×óÁÐÑÓ³Ù¸´ÖÆ
                        C[beta + index_2 - 1][2 * l_index + 1] = C[beta - 1][2 * l_index + 1];
                    }
                }
                index_1 = lambda_offset[layer];
                index_2 = lambda_offset[layer + 1];
                for (int beta = index_1; beta < 2 * index_1; ++beta) {
                    C[beta + index_1 - 1][2 * l_index] = (C[beta - 1][2 * lazy_copy[layer][l_index]] + C[beta - 1][2 * l_index + 1]) % 2; // ×óÁÐÑÓ³Ù¸´ÖÆ
                    C[beta + index_2 - 1][2 * l_index] = C[beta - 1][2 * l_index + 1];
                }

            }

        }

        //lazy_copy----------------
        if (phi < N - 1) {
            for (int i_layer = 0; i_layer <= llr_layer_vec[phi + 1]; ++i_layer) {
                for (int l_index = 0; l_index < L; ++l_index) {
                    lazy_copy[i_layer][l_index] = l_index;
                }
            }
        }
    }

    //crc selection
    std::vector<int> path_ordered(L);
    for (int i = 0; i < L; ++i) {
        path_ordered[i] = i; // Initialize with path indices
    }
    std::sort(path_ordered.begin(), path_ordered.end(),
        [this](int a, int b) { return PM[a] < PM[b]; });
    _U = checkvalue;
    for (int l_index = 0; l_index < L; ++l_index) {
        int path_num = path_ordered[l_index];
        if (crc32 == crc->calculate(u[path_num])) {
            _U = u[path_num];
            IsSucessful = true;
            break;
        }
    }
    return IsSucessful;
}

template <typename B, typename R>
void DecoderPolar<B, R>
::get_lambdaoffset(){
    int numBits = static_cast<int>(std::log2(this->N)) + 1;
    int powerOfTwo = 1;
    for (int i = 0; i < numBits; ++i) {
        this->lambda_offset.push_back(powerOfTwo);
        powerOfTwo *= 2;
    }
}

template <typename B, typename R>
void DecoderPolar<B, R>
::get_bit_layer(){
    for (int phi = 0; phi < this->N; ++phi) {
        int psi = phi / 2;
        int layer = 0;
        while (psi % 2 == 1) {
            psi = psi / 2;
            layer++;
        }
        this->bit_layer_vec[phi] = layer;
    }
}

template <typename B, typename R>
void DecoderPolar<B, R>
::get_llr_layer(){
    for (int phi = 1; phi < this->N; ++phi) {
        int psi = phi;
        int layer = 0;
        while (psi % 2 == 0) {
            psi = psi / 2;
            layer++;
        }
        this->llr_layer_vec[phi] = layer;
    }
}

template<typename B, typename R>
R DecoderPolar<B, R>
::f(R L1, R L2){
    return std::copysign(1.0, L1) * std::copysign(1.0, L2) * std::min(std::fabs(L1), std::fabs(L2));
}

template<typename B, typename R>
void DecoderPolar<B, R>
::bittollr(const std::vector<B>& U, std::vector<R>& llr){
    for (unsigned i = 0; i < llr.size(); i++) {
        llr[i] = U[i] ? llr1 : llr0;
    }
}

template<typename B, typename R>
void DecoderPolar<B, R>
::frozenbitoperation(int phi, int phi_mod_2, const std::vector<B>& checkvalue){
    for (int l_index = 0; l_index < L; ++l_index) {
        if (activepath[l_index] == 0)
            continue;

        u[l_index][phi] = checkvalue[phi];
        if (P[0][l_index] < 0 && checkvalue[phi] == 0) {
            PM[l_index] = PM[l_index] - P[0][l_index];
        }
        if (P[0][l_index] > 0 && checkvalue[phi] == 1) {
            PM[l_index] = PM[l_index] + P[0][l_index];
        }
        if (phi_mod_2 == 0) {
            C[0][2 * l_index] = checkvalue[phi];
        }
        else {
            C[0][2 * l_index + 1] = checkvalue[phi];
        }
    }
}

template<typename B, typename R>
void DecoderPolar<B, R>
::unfrozenbitoperation(int phi, int phi_mod_2) {
    std::fill(PM_pair.begin(), PM_pair.end(), std::numeric_limits<R>::max());
    for (int l_index = 0; l_index < L; ++l_index) {
        if (activepath[l_index] == 0) {
            continue;
        }
        if (P[0][l_index] >= 0) {
            PM_pair[l_index] = PM[l_index];
            PM_pair[l_index + L] = PM[l_index] + P[0][l_index];
        }
        else {
            PM_pair[l_index] = PM[l_index] - P[0][l_index];
            PM_pair[l_index + L] = PM[l_index];
        }
    }
    int middle = std::min(2 * std::accumulate(activepath.begin(), activepath.end(), 0), L);
    std::iota(PM_inx.begin(), PM_inx.end(), 0);
    std::sort(PM_inx.begin(), PM_inx.end(), [&](size_t a, size_t b) {
        return PM_pair[a] < PM_pair[b];
        });
    std::fill(compare.begin(), compare.end(), 0);
    for (size_t i = 0; i < middle; ++i) {
        compare[PM_inx[i]] = 1;
    }
    std::stack<int> kill_index;
    for (int i = 0; i < L; ++i) {
        if (compare[i] == 0 && compare[i + L] == 0) {
            activepath[i] = 0;
            kill_index.push(i);
        }
    }
    for (int l_index = 0; l_index < L; ++l_index) {
        if (activepath[l_index] == 0) {
            continue;
        }
        int path_state = compare[l_index] * 2 + compare[l_index + L];
        switch (path_state) {
        case 1:
            u[l_index][phi] = 1;
            C[0][2 * l_index + phi_mod_2] = 1;
            PM[l_index] = PM_pair[l_index + L];
            break;
        case 2:
            u[l_index][phi] = 0;
            C[0][2 * l_index + phi_mod_2] = 0;
            PM[l_index] = PM_pair[l_index];
            break;
        case 3: {
            int index = kill_index.top();
            kill_index.pop();
            activepath[index] = 1;
            // Lazy copy
            for (int i = 0; i < lazy_copy.size(); ++i) {
                lazy_copy[i][index] = lazy_copy[i][l_index];
            }
            std::copy(u[l_index].begin(), u[l_index].begin() + phi, u[index].begin());
            u[l_index][phi] = 0;
            u[index][phi] = 1;
            C[0][2 * l_index + phi_mod_2] = 0;
            C[0][2 * index + phi_mod_2] = 1;
            PM[l_index] = PM_pair[l_index];
            PM[index] = PM_pair[l_index + L];
            break;
        }
        }
    }
}


template <typename B = int, typename R = float>
class DecoderPcPolar : public DecoderPolar<B> {
private:

public:
    DecoderPcPolar(const int& N, const int& L, const float QBER, const std::shared_ptr<CRC<B>>& crc, const std::shared_ptr<FrozenSiteGenerator>& fbg,int p=9);

    ~DecoderPcPolar();
protected:
    int p;
    std::vector<std::vector<B>>   reg;
    int  idx_reg;
    virtual void decodeinit() override;
    virtual void unfrozenbitoperation(int phi, int phi_mod_2)override;
    virtual void frozenbitoperation(int phi, int phi_mod_2, const std::vector<B>& checkvalue)override;
};

template<typename B, typename R>
DecoderPcPolar<B, R>
::DecoderPcPolar(const int& N, const int& L, const float QBER, const std::shared_ptr<CRC<B>>& crc, const std::shared_ptr<FrozenSiteGenerator>& fbg, int p)
    :DecoderPolar<B,R>(N, L, QBER, crc, fbg), p(p),reg(L, std::vector<B>(p, 0)), idx_reg(0){
}

template<typename B, typename R>
DecoderPcPolar<B, R>
::~DecoderPcPolar() {

}

template<typename B, typename R>
void DecoderPcPolar<B, R>
::frozenbitoperation(int phi, int phi_mod_2, const std::vector<B>& checkvalue) {
    for (int l_index = 0; l_index < this->L; ++l_index) {
        if (this->activepath[l_index] == 0)
            continue;

        if (checkvalue[phi] == this->reg[l_index][idx_reg]) {
            this->PM[l_index] = this->PM[l_index] - (this->P[0][l_index] < 0 ? this->P[0][l_index] : 0);
            this->u[l_index][phi] = 0;
            this->C[0][2 * l_index + phi_mod_2] = 0;
        }
        else {
            this->PM[l_index] = this->PM[l_index] + (this->P[0][l_index] > 0 ? this->P[0][l_index] : 0);
            this->u[l_index][phi] = 1;
            this->C[0][2 * l_index + phi_mod_2] = 1;
            //this->reg[l_index][idx_reg] = 1- this->reg[l_index][idx_reg] ;
        }
    }
    this->idx_reg++;
    if (this->idx_reg == p) { this->idx_reg = 0; }
}

template<typename B, typename R>
void DecoderPcPolar<B, R>
::decodeinit(){
    idx_reg = 0; 
    std::fill(reg.begin(), reg.end(), std::vector<B>(p, 0));
}

template<typename B, typename R>
void DecoderPcPolar<B, R>
::unfrozenbitoperation(int phi, int phi_mod_2) {
    std::fill(this->PM_pair.begin(), this->PM_pair.end(), std::numeric_limits<R>::max());
    for (int l_index = 0; l_index < this->L; ++l_index) {
        if (this->activepath[l_index] == 0) {
            continue;
        }
        if (this->P[0][l_index] >= 0) {
            this->PM_pair[l_index] = this->PM[l_index];
            this->PM_pair[l_index + this->L] = this->PM[l_index] + this->P[0][l_index];
        }
        else {
            this->PM_pair[l_index] = this->PM[l_index] - this->P[0][l_index];
            this->PM_pair[l_index + this->L] = this->PM[l_index];
        }
    }
    int middle = std::min(2 * std::accumulate(this->activepath.begin(), this->activepath.end(), 0), this->L);
    std::iota(this->PM_inx.begin(), this->PM_inx.end(), 0);
    std::sort(this->PM_inx.begin(), this->PM_inx.end(), [&](size_t a, size_t b) {
        return this->PM_pair[a] < this->PM_pair[b];
        });
    std::fill(this->compare.begin(), this->compare.end(), 0);
    for (size_t i = 0; i < middle; ++i) {
        this->compare[this->PM_inx[i]] = 1;
    }
    std::stack<int> kill_index;
    for (int i = 0; i < this->L; ++i) {
        if (this->compare[i] == 0 && this->compare[i + this->L] == 0) {
            this->activepath[i] = 0;
            kill_index.push(i);
        }
    }
    for (int l_index = 0; l_index < this->L; ++l_index) {
        if (this->activepath[l_index] == 0) {
            continue;
        }
        int path_state = this->compare[l_index] * 2 + this->compare[l_index + this->L];
        switch (path_state) {
            case 1:
                this->u[l_index][phi] = 1;
                this->C[0][2 * l_index + phi_mod_2] = 1;
                this->PM[l_index] = this->PM_pair[l_index + this->L];
                this->reg[l_index][idx_reg] = 1 - this->reg[l_index][idx_reg];
                break;
            case 2:
                this->u[l_index][phi] = 0;
                this->C[0][2 * l_index + phi_mod_2] = 0;
                this->PM[l_index] = this->PM_pair[l_index];
                break;
            case 3: {
                int index = kill_index.top();
                kill_index.pop();
                this->activepath[index] = 1;
                // Lazy copy
                for (int i = 0; i < this->lazy_copy.size(); ++i) {
                    this->lazy_copy[i][index] = this->lazy_copy[i][l_index];
                }
                std::copy(this->u[l_index].begin(), this->u[l_index].begin() + phi, this->u[index].begin());
                this->u[l_index][phi] = 0;
                this->u[index][phi] = 1;
                this->C[0][2 * l_index + phi_mod_2] = 0;
                this->C[0][2 * index + phi_mod_2] = 1;
                this->PM[l_index] = this->PM_pair[l_index];
                this->PM[index] = this->PM_pair[l_index + this->L];
                this->reg[index] = this->reg[l_index];
                this->reg[index][idx_reg] = 1 - this->reg[index][idx_reg];
                break;
            }
        }
    }
    this->idx_reg++;
    if (this->idx_reg == p) { this->idx_reg = 0; }
}


template <typename B = int, typename R = float>
class DecoderSegPcPolar : public DecoderPcPolar<B> {
private:

public:
    DecoderSegPcPolar(const int& N, const int& L, const float QBER, const std::shared_ptr<CRC<B>>& crc, const std::shared_ptr<FrozenSiteGenerator>& fbg);

    ~DecoderSegPcPolar();
protected:
    virtual void unfrozenbitoperation(int phi, int phi_mod_2)override;
    virtual void decodeinit() override;
    bool iserror;
};

template<typename B, typename R>
DecoderSegPcPolar<B, R>
::DecoderSegPcPolar(const int& N, const int& L, const float QBER, const std::shared_ptr<CRC<B>>& crc, const std::shared_ptr<FrozenSiteGenerator>& fbg)
:DecoderPcPolar<B, R>(N,L,QBER,crc,fbg) , iserror(false){

}

template<typename B, typename R>
DecoderSegPcPolar<B, R>
::~DecoderSegPcPolar(){

}

template<typename B, typename R>
void DecoderSegPcPolar<B, R>
::unfrozenbitoperation(int phi, int phi_mod_2){
    if (iserror == false) {
        std::fill(this->PM_pair.begin(), this->PM_pair.end(), std::numeric_limits<R>::max());
        for (int l_index = 0; l_index < this->L; ++l_index) {
            if (this->activepath[l_index] == 0) {
                continue;
            }
            if (this->P[0][l_index] >= 0) {
                this->PM_pair[l_index] = this->PM[l_index];
                this->PM_pair[l_index + this->L] = this->PM[l_index] + this->P[0][l_index];
            }
            else {
                this->PM_pair[l_index] = this->PM[l_index] - this->P[0][l_index];
                this->PM_pair[l_index + this->L] = this->PM[l_index];
            }
        }
        int middle = std::min(2 * std::accumulate(this->activepath.begin(), this->activepath.end(), 0), this->L);
        std::iota(this->PM_inx.begin(), this->PM_inx.end(), 0);
        std::sort(this->PM_inx.begin(), this->PM_inx.end(), [&](size_t a, size_t b) {
            return this->PM_pair[a] < this->PM_pair[b];
            });
        std::fill(this->compare.begin(), this->compare.end(), 0);
        for (size_t i = 0; i < middle; ++i) {
            this->compare[this->PM_inx[i]] = 1;
        }
        std::stack<int> kill_index;
        for (int i = 0; i < this->L; ++i) {
            if (this->compare[i] == 0 && this->compare[i + this->L] == 0) {
                this->activepath[i] = 0;
                kill_index.push(i);
            }
        }

        if (this->fbg->errorsegsite < 3) {
            if (phi == ((this->fbg)->SegSite)[this->fbg->errorsegsite]) {
                for (int l_index = 0; l_index < this->L; ++l_index) {
                    if (this->activepath[l_index] == 0) {
                        continue;
                    }
                    this->u[l_index][phi] = 0;
                    if (this->crc->sitecrccalculate(this->u[l_index], phi) == this->fbg->crc_seg[this->fbg->errorsegsite]) {
                        this->C[0][2 * l_index + phi_mod_2] = 0;
                        this->PM[l_index] = this->PM_pair[l_index];
                        std::fill(this->activepath.begin(), this->activepath.end(), 0);
                        this->activepath[l_index] = 1;
                        this->fbg->errorsegsite++;
                        break;
                    }
                    this->u[l_index][phi] = 1;
                    if (this->crc->sitecrccalculate(this->u[l_index], phi) == this->fbg->crc_seg[this->fbg->errorsegsite]) {
                        this->C[0][2 * l_index + phi_mod_2] = 1;
                        this->PM[l_index] = this->PM_pair[l_index + this->L];
                        this->reg[l_index][this->idx_reg] = 1 - this->reg[l_index][this->idx_reg];
                        std::fill(this->activepath.begin(), this->activepath.end(), 0);
                        this->activepath[l_index] = 1;
                        this->fbg->errorsegsite++;
                        break;
                    }
                    this->activepath[l_index] = 0;
                }
                if (std::accumulate(this->activepath.begin(), this->activepath.end(), 0) == 0) {
                    iserror = true;
                }
            }
            else {
                for (int l_index = 0; l_index < this->L; ++l_index) {
                    if (this->activepath[l_index] == 0) {
                        continue;
                    }
                    int path_state = this->compare[l_index] * 2 + this->compare[l_index + this->L];
                    switch (path_state) {
                    case 1:
                        this->u[l_index][phi] = 1;
                        this->C[0][2 * l_index + phi_mod_2] = 1;
                        this->PM[l_index] = this->PM_pair[l_index + this->L];
                        this->reg[l_index][this->idx_reg] = 1 - this->reg[l_index][this->idx_reg];
                        break;
                    case 2:
                        this->u[l_index][phi] = 0;
                        this->C[0][2 * l_index + phi_mod_2] = 0;
                        this->PM[l_index] = this->PM_pair[l_index];
                        break;
                    case 3: {
                        int index = kill_index.top();
                        kill_index.pop();
                        this->activepath[index] = 1;
                        // Lazy copy
                        for (int i = 0; i < this->lazy_copy.size(); ++i) {
                            this->lazy_copy[i][index] = this->lazy_copy[i][l_index];
                        }
                        std::copy(this->u[l_index].begin(), this->u[l_index].begin() + phi, this->u[index].begin());
                        this->u[l_index][phi] = 0;
                        this->u[index][phi] = 1;
                        this->C[0][2 * l_index + phi_mod_2] = 0;
                        this->C[0][2 * index + phi_mod_2] = 1;
                        this->PM[l_index] = this->PM_pair[l_index];
                        this->PM[index] = this->PM_pair[l_index + this->L];
                        this->reg[index] = this->reg[l_index];
                        this->reg[index][this->idx_reg] = 1 - this->reg[index][this->idx_reg];
                        break;
                    }
                    }
                }
            }
        }
        else {
            for (int l_index = 0; l_index < this->L; ++l_index) {
                if (this->activepath[l_index] == 0) {
                    continue;
                }
                int path_state = this->compare[l_index] * 2 + this->compare[l_index + this->L];
                switch (path_state) {
                case 1:
                    this->u[l_index][phi] = 1;
                    this->C[0][2 * l_index + phi_mod_2] = 1;
                    this->PM[l_index] = this->PM_pair[l_index + this->L];
                    this->reg[l_index][this->idx_reg] = 1 - this->reg[l_index][this->idx_reg];
                    break;
                case 2:
                    this->u[l_index][phi] = 0;
                    this->C[0][2 * l_index + phi_mod_2] = 0;
                    this->PM[l_index] = this->PM_pair[l_index];
                    break;
                case 3: {
                    int index = kill_index.top();
                    kill_index.pop();
                    this->activepath[index] = 1;
                    // Lazy copy
                    for (int i = 0; i < this->lazy_copy.size(); ++i) {
                        this->lazy_copy[i][index] = this->lazy_copy[i][l_index];
                    }
                    std::copy(this->u[l_index].begin(), this->u[l_index].begin() + phi, this->u[index].begin());
                    this->u[l_index][phi] = 0;
                    this->u[index][phi] = 1;
                    this->C[0][2 * l_index + phi_mod_2] = 0;
                    this->C[0][2 * index + phi_mod_2] = 1;
                    this->PM[l_index] = this->PM_pair[l_index];
                    this->PM[index] = this->PM_pair[l_index + this->L];
                    this->reg[index] = this->reg[l_index];
                    this->reg[index][this->idx_reg] = 1 - this->reg[index][this->idx_reg];
                    break;
                }
                }
            }
        }

        
        this->idx_reg++;
        if (this->idx_reg == this->p) { this->idx_reg = 0; }
    }
}


template<typename B, typename R>
void DecoderSegPcPolar<B, R>
::decodeinit(){
    this->idx_reg = 0;
    std::fill(this->reg.begin(), this->reg.end(), std::vector<B>(this->p, 0));
    this->iserror = false;
    this->fbg->errorsegsite = 0;
}
#endif // DECODER_HPP_