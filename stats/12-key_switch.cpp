#include "../rlwe.h"
#include "../rgsw.h"
#include "../ksw_key.h"
#include "../operations.h"
#include "../param.h"

#define NB_TESTS 25

typedef CirculantRing<Zt, P1> MyRpt;
int main()
{
    /*
    double cumul = 0;
    for (size_t i = 0; i < 1000 ; ++i)
    {
        int64_t value 
        cumul += (int64_t)((Rz::gaussian_sample(8.0)).get_data()[0]);
    }
    std::cout << cumul / 1000 << std::endl;
    */
    int seed = time(NULL);
    srand(seed);
    //std::cout << "seed " << seed << std::endl;
    print_param();

    RLWE<Rz, P1> in;
    double noise = 0;
    double cumul_noise = 0;

    Rz *s = new Rz[P1];
    for (size_t i = 0 ; i < P1 ; ++i)
        s[i] = Rz::sample_s(DENSITY_KEY);
    Rz *s2 = new Rz[N];
    for (size_t i = 0 ; i < N ; ++i)
        s2[i] = Rz::sample_s(DENSITY_KEY_SMALL);
    for (size_t idx = 0 ; idx < NB_TESTS ; ++idx)
    {
        KSWKeyLWE *S_lwe = new KSWKeyLWE(s, s2, VARIANCE_KSWLWE);
        
        CirculantRing<Zt, 1, 1> m = CirculantRing<Zt, 1, 1>::uniform_sample();
        in.encrypt(s, m, Qp, T, std::pow(2, 91.22));
        RLWE<Rz, N> out = in.key_switch(*S_lwe);

        noise = out.noise(s2, Qp, T);
        cumul_noise += noise;
        delete S_lwe;
    }
    std::cout << cumul_noise << " " << cumul_noise / NB_TESTS << " = 2^" << std::log2(cumul_noise / NB_TESTS) << std::endl;
    delete[] s;
    delete[] s2;

    return 0;
}
