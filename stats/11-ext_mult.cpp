
#include "../rlwe.h"
#include "../rgsw.h"
#include "../ksw_key.h"
#include "../operations.h"
#include "../param.h"

#define NB_TESTS 1000

typedef CirculantRing<Zt, P1> MyRpt;
int main()
{
    int seed = time(NULL);
    srand(seed);
    //std::cout << "seed " << seed << std::endl;
    print_param();

    Rp1GSW::G_init();
    double noise = 0;
    double cumul_noise = 0;

    for (size_t idx = 0 ; idx < NB_TESTS ; ++idx)
    {
        Rp1 s[1];
        s[0] = Rp1::sample_s(DENSITY_KEY);
        /** Testing basic encryption and decryption **/
        size_t m = rand() % (P1);
        Zt Xm_coefs[P1];
        for (size_t j = 0 ; j < P1 ; ++j) Xm_coefs[j] = 0;
        Xm_coefs[m] = 1;
        MyRpt Xm(Xm_coefs);
        //std::cout << "Xm\t" << Xm << std::endl;
        Rp1LWE rlwe;
        rlwe.encrypt(s, Xm, Q, T, 1);

        size_t mp = rand() % (P1);
        Zt Xmp_coefs[P1];
        for (size_t j = 0 ; j < P1 ; ++j) Xmp_coefs[j] = 0;
        Xmp_coefs[mp] = 1;
        MyRpt Xmp(Xmp_coefs);
        //std::cout << "Xmp\t" << Xmp << std::endl;
        Rp1GSW rgsw;
        rgsw.encrypt(s, Xmp, VARIANCE_ACC);

        /** Testing ExtMult **/
        RLWE<Rp1, 1> rlwe2 = rlwe.ext_mult(rgsw);
        Zt Xm_mp_coefs[P1];
        for (size_t j = 0 ; j < P1 ; ++j) Xm_mp_coefs[j] = 0;
        Xm_mp_coefs[(m + mp) % P1] = 1;
        MyRpt Xm_mp(Xm_mp_coefs);
        //std::cout << "Noise after: " <<  rlwe2.noise(s, Q, T) << std::endl;
        noise = rlwe2.noise(s, Q, T);
        cumul_noise += noise;
    }
    std::cout << cumul_noise << " " << cumul_noise / NB_TESTS / P1 << " = 2^" << std::log2(cumul_noise / NB_TESTS / P1) << std::endl;

    return 0;
}
