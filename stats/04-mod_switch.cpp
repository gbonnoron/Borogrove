#include "../operations.h"
#include "../rlwe.h"
#include "../param.h"

#define NB_TESTS 1000

int main()
{
    int seed = time(NULL);
    srand(seed);
    //std::cout << "seed " << seed << std::endl;
    print_param();

    Rp1LWE rlwe;
    RLWE<Rp1_crt, 1> rlwe_ms;
    CirculantRing<Zt, P1> m;
    Rp1 s[P1];
    Rp1_crt s_ms[P1];

    double noise = 0;
    double cumul_noise = 0;

    for (size_t idx = 0 ; idx < NB_TESTS ; ++idx)
    {
        m = CirculantRing<Zt, P1>::uniform_sample();
        for (size_t i = 0 ; i < P1 ; ++i)
        {
            s[i] = Rp1::sample_s(DENSITY_KEY);
            s_ms[i] = s[i];
        }
        rlwe.encrypt(s, m, Q, T, std::pow(2, 38.50));
        //rlwe.encrypt(s, m, Q, T, 0);
        //rlwe.encrypt(s, m, Q, T, 32);
        //std::cout << "Noise before: " << lwe.noise(s, Q, T) << std::endl;

        rlwe_ms = rlwe.mod_switch<Rp1_crt, Qcrt, Q>();

        //std::cout << "Noise after: " << rlwe_ms.noise(s_ms, Qp, T) << std::endl;
        noise = rlwe_ms.noise(s_ms, Qcrt, T);
        cumul_noise += noise;
    }
    std::cout << cumul_noise << " " << cumul_noise / NB_TESTS / P1 << " = 2^" << std::log2(cumul_noise / NB_TESTS / P1) << std::endl;

    return 0;
}
