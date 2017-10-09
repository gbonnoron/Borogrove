#include "../rlwe.h"
#include "../rgsw.h"
#include "../ksw_key.h"
#include "../operations.h"
#include "../param.h"

#define NB_TESTS 100

int main()
{
    int seed = time(NULL);
    srand(seed);
    //std::cout << "seed " << seed << std::endl;
    print_param();

    Rp1GSW::G_init();

    Rp1 s[1];
    unsigned int l = N;
    KSWKeyRp1 S[P1];
    Rp1GSW C[l+1];
    Zqp x[l];
    Rz s_x[l];
    Zp1 x2[l];
    Zp1 y[l];

    double noise = 0;
    double cumul_noise = 0;

    s[0] = Rp1::sample_s(DENSITY_KEY);
    gen_keyswitching_keys<Rp1, P1>(S, s[0], VARIANCE_ACC);
    for (unsigned int i = 0 ; i < l ; ++i)
    {
        x[i] = rand();
        x2[i] = x[i] % P1;
        s_x[i] = x[i];
    }
    gen_bootstrapping_keys<Rp1, N, P1>(Q, T, C, s_x, s[0], VARIANCE_ACC);
    for (size_t idx = 0 ; idx < NB_TESTS ; ++idx)
    {
        for (unsigned int i = 0 ; i < l ; ++i)
            y[i] = rand();

        //std::cout << "Noise before: " << C[0].noise(s, Q, T) << std::endl;
        Rp1LWE c = ext_exp_inner<Rp1, Zp1, P1>(Q, T, l, y, C, S);

        //std::cout << "Noise after: " <<  c.noise(s, Q, T) << std::endl;
        noise = c.noise(s, Q, T);
        cumul_noise += noise;
    }
    std::cout << cumul_noise << " " << cumul_noise / NB_TESTS / P1 << " = 2^" << std::log2(cumul_noise / NB_TESTS / P1) << std::endl;

    return 0;
}
