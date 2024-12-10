#ifndef ENCODER_HPP_
#define ENCODER_HPP_

#include <vector>

template <typename B = int>
class Encoder {
public:
    Encoder() {}

    virtual ~Encoder() = default;

    virtual void encode(const std::vector<B>& Ksa, std::vector<B>& U) = 0;
};

template <typename B = int>
class EncoderPolar : public Encoder<B>
{
public:
    EncoderPolar() {}

    ~EncoderPolar() {}

    void encode(const std::vector<B>& Ksa, std::vector<B>& U) override;
};

template<typename B>
void EncoderPolar<B>
::encode(const std::vector<B>& Ksa, std::vector<B>& U) {
    int N = Ksa.size();
    U = Ksa;
    for (auto k = (N >> 1); k > 0; k >>= 1)
        for (auto j = 0; j < N; j += 2 * k)
            for (auto i = 0; i < k; i++)
                U[j + i] = U[j + i] ^ U[k + j + i];
}

#endif // ENCODER_HPP_