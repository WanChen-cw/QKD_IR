#ifndef FROZENSITEGENERATOR_HPP_
#define FROZENSITEGENERATOR_HPP_

#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <algorithm>
#include <bitset>

int count0 = 0;
int count1 = 0;
int count2 = 0;
int count3 = 0;

class FrozenSiteGenerator
{
private:
    int N; /*!< Codeword size (or frame size). */
    int Nf;
    std::string path;
    std::vector<uint32_t> best_channels; /*!< The best channels in a codeword sorted by descending order. */
public:
    std::vector<bool> frozen_bits;

    FrozenSiteGenerator(int N, const std::string& path = " ");

    ~FrozenSiteGenerator();

    void generate(int Nf);//init frozenbit site

    void addfb(int t);//add frozenbit site

    int addfb_s(int t);//add frozenbit site

    void SegSiteGenerate(int Ns, int Nf_min, int Nf_max);

    std::vector<int> SegSite;
    std::vector<uint32_t> crc_seg;
    int errorsegsite;

protected:
    void get_best_channels();
};

FrozenSiteGenerator
::FrozenSiteGenerator(int N, const std::string& path)
    :N(N),Nf(0), path(path), frozen_bits(N,true), errorsegsite(0){
    this->get_best_channels();
}

FrozenSiteGenerator
::~FrozenSiteGenerator(){
}

void FrozenSiteGenerator
::generate(int nf){
    if (Nf != nf) {
        std::fill(frozen_bits.begin(), frozen_bits.end(), false);
        for (auto i = N-1; i > N-1-nf; i--){
            frozen_bits[best_channels[i]] = true;
        }
        Nf = nf;
    }
}

void FrozenSiteGenerator
::addfb(int t) {
    for (auto i = N-Nf-1; i > N-Nf-1-t; i--) {
        frozen_bits[best_channels[i]] = true;
    }
    Nf += t;
}

int FrozenSiteGenerator
::addfb_s(int t){
    int temp = 0;
    for (auto i = N - Nf - 1; i > N - Nf - 1 - t; i--) {
        if (errorsegsite == 0) {
            frozen_bits[best_channels[i]] = true;
            temp++;
        }
        else if (best_channels[i]> SegSite[errorsegsite-1]) {
            frozen_bits[best_channels[i]] = true;
            temp++;
        }  
    }
    Nf += t;
    switch (errorsegsite)
    {
    case 0:     
        count0++;
        break;
    case 1:
        count1++;
        break;
    case 2:
        count2++;
        break;
    case 3:
        count3++;
        break;
    default:
        break;
    }
    return temp;
}

void FrozenSiteGenerator
::SegSiteGenerate(int Ns, int Nf_min, int Nf_max){
    int start = N - Nf_max; 
    int end = N - Nf_min;   
    if (start < 0 || end < 0 || start > end || end >= best_channels.size()) {
        std::cerr << "error:segsitegenerete\n" << std::endl;
    }

    std::vector<uint32_t> segment(best_channels.begin() + start, best_channels.begin() + end + 1);
    std::sort(segment.begin(), segment.end());

    int segment_size = segment.size();
   
    
    /////////////////1：1：1：1 （8：4：2：1）      fn: 791 | Nf+Ncrc*4: 10659 | f: 1.1499 | rounds: 3.23262 | time: 55 s5 s
    //int step = segment_size / Ns;  // 每段的大小（整数部分）
    //std::vector<int> boundaries(Ns-1);
    //for (int i = 0; i < Ns-1; ++i) {
    //    int end_idx = step* (i+1);
    //    int idx = (segment[end_idx]+ segment[end_idx+1])/2;
    //    while (idx == segment[end_idx]||frozen_bits[idx]==true) {
    //        end_idx++;
    //        idx = (segment[end_idx] + segment[end_idx + 1]) / 2;
    //    }
    //    boundaries[i]= idx;
    //}

    /////////////////1：2：2：3（4：4：3：2 ）fn: 660 | Nf+Ncrc*4: 10663 | f: 1.15034 | rounds: 3.30152 | time: 45 s s    
    //int step = segment_size / Ns/2;  // 每段的大小（整数部分）
    //std::vector<int> boundaries(Ns - 1);
    //for (int i = 0; i < Ns - 1; ++i) {
    //    int end_idx = step * (i*2 + 1);
    //    int idx = (segment[end_idx] + segment[end_idx + 1]) / 2;
    //    while (idx == segment[end_idx] || frozen_bits[idx] == true) {
    //        end_idx++;
    //        idx = (segment[end_idx] + segment[end_idx + 1]) / 2;
    //    }
    //    boundaries[i] = idx;
    //}

    /////////1；1：2：4 （4：3：3：3）  fn: 1114 | Nf+Ncrc*4: 10645.4 | f: 1.14844 | rounds: 3.22083 | time: 85 s
    int step = segment_size / 8;  // 每段的大小（整数部分）
    std::vector<int> boundaries(Ns - 1);

        int end_idx = step * 1;
        int idx = (segment[end_idx] + segment[end_idx + 1]) / 2;
        while (idx == segment[end_idx] || frozen_bits[idx] == true) {
            end_idx++;
            idx = (segment[end_idx] + segment[end_idx + 1]) / 2;
        }
        boundaries[0] = idx;

        end_idx = step * 2;
        idx = (segment[end_idx] + segment[end_idx + 1]) / 2;
        while (idx == segment[end_idx] || frozen_bits[idx] == true) {
            end_idx++;
            idx = (segment[end_idx] + segment[end_idx + 1]) / 2;
        }
        boundaries[1] = idx;

        end_idx = step * 4;
        idx = (segment[end_idx] + segment[end_idx + 1]) / 2;
        while (idx == segment[end_idx] || frozen_bits[idx] == true) {
            end_idx++;
            idx = (segment[end_idx] + segment[end_idx + 1]) / 2;
        }
        boundaries[2] = idx;


    SegSite= boundaries;
}

void FrozenSiteGenerator
::get_best_channels(){
    std::ifstream inputFile(this->path);
    if (!inputFile.is_open()) {
        std::cout << "Failed to open channel parameter file: " << this->path << std::endl;
        return;
    }
    uint32_t value;
    while (inputFile >> value) {
        this->best_channels.push_back(value);
    }
    inputFile.close();
}
#endif /* FROZENSITEGENERATOR_HPP_ */
