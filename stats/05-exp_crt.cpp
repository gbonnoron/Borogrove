#include "../rlwe.h"
#include "../operations.h"
#include "../param.h"

#define NB_TESTS 100

int main()
{
    int seed = time(NULL);
    srand(seed);
    //std::cout << "seed " << seed << std::endl;
    print_param();

    size_t m;
    Zt Xmp_coefs[P1], Ymq_coefs[P2];
    Rp1 s_p_tmp[1];
    Rp2 s_q_tmp[1];
    Rp12 s_pq_tmp[3];
    Rp1_crt s_p[1];
    Rp2_crt s_q[1];
    Rp12_crt s_pq[3];
    RLWE<Rp1_crt, 1> c_p;
    RLWE<Rp2_crt, 1> c_q;
    RLWE<Rp12_crt, 3> c_pq;

    double noise = 0;
    double cumul_noise = 0;

    for (size_t idx = 0 ; idx < NB_TESTS ; ++idx)
    {
        m = rand() % (P1*P2);

        for (size_t j = 0 ; j < P1 ; ++j) Xmp_coefs[j] = 0;
        Xmp_coefs[m % P1] = 1;
        CirculantRing<Zt, P1> Xmp(Xmp_coefs);

        for (size_t j = 0 ; j < P2 ; ++j) Ymq_coefs[j] = 0;
        Ymq_coefs[m % P2] = 1;
        CirculantRing<Zt, P2> Ymq(Ymq_coefs);

        s_p[0] = Rp1_crt::sample_s(DENSITY_KEY);
        s_q[0] = Rp2_crt::sample_s(DENSITY_KEY);
        s_p_tmp[0] = s_p[0];
        s_q_tmp[0] = s_q[0];

        c_p.encrypt(s_p, Xmp, Qcrt, T, std::pow(2, 6.41));
        c_q.encrypt(s_q, Ymq, Qcrt, T, std::pow(2, 6.41));
        //std::cout << "Noise before: " << c_p.noise(s_p, Qp, T) << std::endl;
        //std::cout << "Noise before: " << c_p.noise(s_q, Qp, T) << std::endl;

        exp_crt(c_pq, c_p, c_q);
        crt_key(s_pq_tmp, s_p_tmp[0], s_q_tmp[0]);
        s_pq[0] = s_pq_tmp[0];
        s_pq[1] = s_pq_tmp[1];
        s_pq[2] = s_pq_tmp[2];

        //std::cout << "Noise after: " << c_pq.noise(s_pq_crt, Qcrt, T) << std::endl;
        //CirculantRing<Zt, P1*P2, FFT_DIM2> dec;
        //c_pq.decrypt(s_pq, dec, Qcrt, T);
        //std::cout << dec << std::endl;
        noise = c_pq.noise(s_pq, Qcrt, T);
        cumul_noise += noise;
    }
    std::cout << cumul_noise << " " << cumul_noise / NB_TESTS / (P1*P2) << " = 2^" << std::log2(cumul_noise / NB_TESTS / (P1*P2)) << std::endl;

    return 0;
}
