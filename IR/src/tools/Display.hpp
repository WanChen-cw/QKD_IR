#include <iostream>
#include <iomanip>
#include <mutex>
#include <chrono> 
#include <cmath>

class Reporter {
private:
    size_t fn;
    size_t fe;
    size_t fnl;
    size_t fel;
    double rounds;
    double addNf_s;
    double f;
    int num_threads;
    bool runcondition;
    mutable std::mutex mtx;
    std::chrono::steady_clock::time_point startTime;

    //parament
    int N;
    int L;
    int t;
    int seed;
    float qber;
    std::string fbgpath;
    int Nf;
    uint32_t crcPolynomial;
    double C;

    std::map<int, int, std::greater<int>> map_t_round;


public:
    Reporter(int N, int L, int seed, float qber, const std::string& fbgpath, int Nf, uint32_t crcPolynomial,int t)
        : fe(0), fn(0), num_threads(0), fnl(100000), fel(100), f(0), rounds(0),  addNf_s(0),runcondition(true),
        N(N), L(L), seed(seed), qber(qber), fbgpath(fbgpath), Nf(Nf), crcPolynomial(crcPolynomial) ,t(t){
        startTime = std::chrono::steady_clock::now();
        C = N * (-qber * log2(qber) - (1 - qber) * log2(1 - qber));
    }

    void set_num_threads(int num_threads) {
        this->num_threads = num_threads;
    }

    void FEAdd(size_t num) {
        std::lock_guard<std::mutex> lock(mtx); 
        fe++;
        fn += num;
        float fer = (fn == 0) ? 0 : (static_cast<float>(fe) / fn) * 100;
        auto currentTime = std::chrono::steady_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(currentTime - startTime).count();

        std::cout << "\rfn: " << fn
            << " | Nf: " << Nf+32
            << " | fe: " << fe
            << " | fer: " << std::fixed << std::setprecision(2) << fer << "%"
            << " | time: " << elapsed << " s"
            << std::flush;
        if (fe >= fel* num_threads) {
            runcondition = false;
        }
    }

    void FNAdd(int num) {
        std::lock_guard<std::mutex> lock(mtx);
        fn++;
        map_t_round[num]++;
        rounds = ((fn - 1) * rounds + num) / fn;
        f = (Nf + 32 + rounds * t) / C;
        auto currentTime = std::chrono::steady_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(currentTime - startTime).count();

        std::cout << "\rfn: " << fn
            << " | Nf+Ncrc: " << Nf + 32 + rounds * t
            << " | f: " << std::setprecision(6) << f
            << " | rounds: " << std::setprecision(6) << rounds + 1
            << " | time: " << elapsed << " s"
            << std::flush;
        if (fn >= fnl) {
            int totalDecodes = 0;  // 总的译码次数
            for (const auto& entry : map_t_round) {
                totalDecodes += entry.second;  // 累加译码次数
            }
            int cumulativeFailures = 0;  // 用来累加当前轮及之前轮的译码次数
            std::cout << "\n" << "每轮译码失败概率：\n";
            for (const auto& entry : map_t_round) {
                int extraRounds = entry.first;  // 额外的译码轮次
                int count = entry.second;  // 对应额外轮次的次数
                // 累加当前轮及之前轮的译码次数
                cumulativeFailures += count;
                // 计算该轮的失败概率
                double failureProbability = static_cast<double>(cumulativeFailures) / totalDecodes;
                std::cout << "轮次 " << extraRounds << " 的失败概率为："
                    << failureProbability * 100 << "%" << std::endl;
            }
            std::cout << std::flush;
            runcondition = false;
        }
    }

    void FNsAdd(int num,int addnf) {
        std::lock_guard<std::mutex> lock(mtx);
        fn++;
        addNf_s= ((fn - 1) * addNf_s + addnf) / fn;
        rounds = ((fn - 1) * rounds + num) / fn;
        f = (Nf + 32*4 + addNf_s) / C;
        auto currentTime = std::chrono::steady_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(currentTime - startTime).count();

        std::cout << "\rfn: " << fn
            << " | Nf+Ncrc*4: " << Nf + 32*4 + addNf_s
            << " | f: " << std::setprecision(6) << f
            << " | rounds: " << std::setprecision(6) << rounds + 1
            << " | time: " << elapsed << " s"
            << std::flush;
        if (fn >= fnl) {
            runcondition = false;
        }
    }

    void printparament(){
        std::cout << "\n---------------------------------" << std::endl;
        std::cout << "N:" << N << std::endl;
        std::cout << "QBER: " << qber << std::endl;
        std::cout << "L: " << L << std::endl;
        std::cout << "seed:" << seed << std::endl;
        std::cout << "t:" << t << std::endl;
        std::cout << "fbgpath: " << fbgpath << std::endl;
        std::cout << "crcPolynomial: " << crcPolynomial << std::endl;
        std::cout << "---------------------------------" << std::endl;
    }

    bool iscontinue() const{
        return runcondition;
    }
};
