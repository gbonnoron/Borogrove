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

    RLWE<Rz, N> Em_short;
    RLWE<CirculantRing<Zp12, 1, 1>, N> ab_short;
    CirculantRing<Zt, 1, 1> m;
    Rz s[N];
    CirculantRing<Zp12, 1, 1> s_ms[N];

    double noise = 0;
    double cumul_noise = 0;

    for (size_t idx = 0 ; idx < NB_TESTS ; ++idx)
    {
        m = CirculantRing<Zt, 1, 1>::uniform_sample();
        for (size_t i = 0 ; i < N ; ++i)
        {
            s[i] = Rz::sample_s(DENSITY_KEY);
            s_ms[i] = s[i];
        }
        Em_short.encrypt(s, m, Qp, T, std::pow(2, 92.48));
        //std::cout << "Noise before: " << lwe.noise(s, Q, T) << std::endl;

        ab_short = Em_short.mod_switch<CirculantRing<Zp12, 1, 1>, P1*P2, Qp>();

        //std::cout << "Noise after: " << lwe_ms.noise(s_ms, P1*P2, T) << std::endl;
        noise = ab_short.noise(s_ms, P1*P2, T);
        cumul_noise += noise;
    }
    std::cout << cumul_noise << " " << cumul_noise / NB_TESTS << " = 2^" << std::log2(cumul_noise / NB_TESTS) << std::endl;

    return 0;
}
