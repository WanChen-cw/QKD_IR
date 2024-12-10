#ifndef SIFTKEYGENERATOR_HPP_
#define SIFTKEYGENERATOR_HPP_

#include <random>
#include <iostream>
#include <fstream>
#include <vector>
#include <memory>

template <typename B=int>
class SiftKeyGenerator {
public:
    SiftKeyGenerator() {}

    virtual ~SiftKeyGenerator() = default; 

    virtual void generate(std::vector<B>& ksa, std::vector<B>& ksb) = 0;
};


template <typename B = int>
class SiftKeyGeneratorContext {
private:
    std::unique_ptr<SiftKeyGenerator<B>> strategy;

public:
    void setStrategy(std::unique_ptr<SiftKeyGenerator<B>> newStrategy) {
        strategy = std::move(newStrategy);
    }

    void generate(std::vector<B>& ksa, std::vector<B>& ksb) {
        if (!strategy) {
            throw std::runtime_error("Strategy not set");
        }
        strategy->generate(ksa, ksb);
    }
};


template <typename B = int>
class RandomSiftKeyGenerator : public SiftKeyGenerator<B> {
private:
    int N;
    float ep;
    std::random_device rd;
    std::mt19937 rd_engine; // Mersenne Twister 19937
#ifdef _MSC_VER
    std::uniform_int_distribution<short> uniform_dist1;
    std::uniform_int_distribution<int> uniform_dist2;
#else
    std::uniform_int_distribution<short> uniform_dist1;
    std::uniform_int_distribution<int> uniform_dist2;
#endif

public:
    RandomSiftKeyGenerator(const int N=1024, const float ep = 0.02,const int seed = 0)
        : SiftKeyGenerator<B>(), N(N), ep(ep),rd(), rd_engine(rd() + seed), uniform_dist1(0, 1), uniform_dist2(0, 100000) {}

    ~RandomSiftKeyGenerator() override {}

    void generate(std::vector<B>& ksa, std::vector<B>& ksb) override {
        for (int i = 0; i < N; i++) {
            ksa[i] = uniform_dist1(rd_engine);
            int random_number = uniform_dist2(rd_engine);
            ksb[i] = random_number > 100000 * ep ? ksa[i] : 1- ksa[i];
        }
    }
};


template <typename B = int>
class FileSiftKeyGenerator : public SiftKeyGenerator<B> {
private:
    int N;
    std::string filenameA;
    std::string filenameB;
    std::ifstream fileA;
    std::ifstream fileB;

public:
    FileSiftKeyGenerator( const std::string& filenameA, const std::string& filenameB,const int N=1024)
        : SiftKeyGenerator<B>(), N(N),filenameA(filenameA), filenameB(filenameB), fileA(filenameA), fileB(filenameB) {
        if (!fileA) {
            throw std::runtime_error("Unable to open fileA");
        }
        if (!fileB) {
            throw std::runtime_error("Unable to open fileB");
        }
    }

    ~FileSiftKeyGenerator() {
        if (fileA.is_open()) {
            fileA.close();
        }
        if (fileB.is_open()) {
            fileB.close();
        }
    }

    void generate(std::vector<B>& ksa, std::vector<B>& ksb) override {
        ksa.clear();
        ksb.clear();
        char chA, chB;
        size_t countA = 0;
        size_t countB = 0;

        while (countA < N && fileA.get(chA) && fileB.get(chB)) {
            if ((chA == '0' || chA == '1') && (chB == '0' || chB == '1')) {
                ksa.push_back(chA - '0');
                ksb.push_back(chB - '0');
                ++countA;
                ++countB;
            }
        }

        if (ksa.size() < N || ksb.size() < N) {
            throw std::runtime_error("Not enough data read from one or both files");
        }
    }

    void reset() {
        if (fileA.is_open()) {
            fileA.clear();
            fileA.seekg(0);
        }
        if (fileB.is_open()) {
            fileB.clear();
            fileB.seekg(0);
        }
    }
};

#endif /* SIFTKEYGENERATOR_HPP_ */