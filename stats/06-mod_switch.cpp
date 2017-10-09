#include "../operations.h"
#include "../rlwe.h"
#include "../param.h"

#define NB_TESTS 10

int main()
{
    int seed = time(NULL);
    srand(seed);
    //std::cout << "seed " << seed << std::endl;
    print_param();

    RLWE<Rp12_crt, 3> rlwe;
    Rp12LWE rlwe_ms;
    CirculantRing<Zt, P1*P2, FFT_DIM2> m;
    Rp12_crt s[3];
    Rp12 s_ms[3];

    double noise = 0;
    double cumul_noise = 0;

    for (size_t idx = 0 ; idx < NB_TESTS ; ++idx)
    {
        m = CirculantRing<Zt, P1*P2, FFT_DIM2>::uniform_sample();
        for (size_t i = 0 ; i < 3 ; ++i)
        {
            s[i] = Rp12_crt::sample_s(DENSITY_KEY);
            s_ms[i] = s[i];
        }
        rlwe.encrypt(s, m, Qcrt, T, std::pow(2, 24.82));
        //std::cout << "Noise before: " << lwe.noise(s, Q, T) << std::endl;

        rlwe_ms = rlwe.mod_switch<Rp12, Qp, Qcrt>();

        //std::cout << "Noise after: " << rlwe_ms.noise(s_ms, Qp, T) << std::endl;
        noise = rlwe_ms.noise(s_ms, Qp, T);
        cumul_noise += noise;
    }
    std::cout << cumul_noise << " " << cumul_noise / NB_TESTS / (P1*P2) << " = 2^" << std::log2(cumul_noise / NB_TESTS / (P1*P2)) << std::endl;

    return 0;
}
