#include "../rlwe.h"
#include "../ksw_key.h"
#include "../operations.h"
#include "../param.h"

#define NB_TESTS 10

int main()
{
    int seed = time(NULL);
    srand(seed);
    //std::cout << "seed " << seed << std::endl;
    print_param();

    size_t m;
    Zt *Zm_coefs = new Zt[P1*P2];
    Fun F;
    Rp12 Rf(F);
    fftw_complex *f = align_alloc<fftw_complex>(ALIGNMENT, FFT_DIM2);
    Rf.compute_fft(f);
    Rp12LWE rlwe_fun;
    Rp12 s_fun[3];
    Rz sp_fun[P1];
    KSWKeyRp12 *S = new KSWKeyRp12();

    LWE lwe_fun;

    double noise = 0;
    double cumul_noise = 0;

    for (size_t j = 0 ; j < P1 ; ++j)
        sp_fun[j] = Rz::sample_s(DENSITY_KEY);
    for (size_t j=0 ; j<3 ; ++j)
        s_fun[j] = Rp12::sample_s(DENSITY_KEY);

    for (size_t idx = 0 ; idx < NB_TESTS ; ++idx)
    {
        m = rand() % (P1*P2);
        for (size_t j = 0 ; j < P1*P2 ; ++j) Zm_coefs[j] = 0;
        Zm_coefs[m % (P1*P2)] = 1;
        CirculantRing<Zt, P1*P2, FFT_DIM2> Zm(Zm_coefs);

        rlwe_fun.encrypt(s_fun, Zm, Qp, T, std::pow(2, 60.82));

        gen_funexpextract_key(S, s_fun, sp_fun, VARIANCE_FUNEXPEXTRACT);
        //std::cout << "Noise before: " << rlwe_fun.noise(s_fun, Qp, T) << std::endl;
        lwe_fun = fun_exp_extract(f, rlwe_fun, *S);

        //std::cout << "Noise after: " <<  lwe_fun.noise(sp_fun, Qp, T) << std::endl;
        noise = lwe_fun.noise(sp_fun, Qp, T);
        cumul_noise += noise;
    }
    std::cout << cumul_noise << " " << cumul_noise / NB_TESTS << " = 2^" << std::log2(cumul_noise / NB_TESTS) << std::endl;
    delete S;
    delete[] Zm_coefs;
    free(f);

    return 0;
}
