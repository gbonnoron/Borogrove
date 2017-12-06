#include <iostream>
#include <string>

#include "../operations.h"
#include "../rlwe.h"
#include "../rgsw.h"
#include "../ksw_key.h"
#include "../param.h"

using namespace std;

double ext_exp_inner(const double var_in, const size_t nb_tests)
{
    Rp1 s[1];
    KSWKeyRp1 S[P1];
    Rp1GSW C[N+1];
    Zqp x[N];
    Rz s_x[N];
    Zp1 x2[N];
    Zp1 y[N];

    double noise = 0;
    double cumul_noise = 0;

    s[0] = Rp1::sample_s(DENSITY_KEY);
    gen_keyswitching_keys<Rp1, P1>(S, s[0], VARIANCE_ACC);
    for (unsigned int i = 0 ; i < N ; ++i)
    {
        x[i] = rand();
        x2[i] = x[i] % P1;
        s_x[i] = x[i];
    }
    gen_bootstrapping_keys<Rp1, N, P1>(Q, T, C, s_x, s[0], VARIANCE_ACC);
    for (size_t idx = 0 ; idx < nb_tests ; ++idx)
    {
        for (unsigned int i = 0 ; i < N ; ++i)
            y[i] = rand();

        //std::cout << "Noise before: " << C[0].noise(s, Q, T) << std::endl;
        Rp1LWE c = ext_exp_inner<Rp1, Zp1, P1>(Q, T, N, y, C, S);

        //std::cout << "Noise after: " <<  c.noise(s, Q, T) << std::endl;
        noise = c.noise(s, Q, T);
        cumul_noise += noise;
    }
    double var_out = log2(cumul_noise / nb_tests / P1);
    return var_out;
}

double mod_switch_acc(const double var_in, const size_t nb_tests)
{
    Rp1LWE rlwe;
    RLWE<Rp1_crt, 1> rlwe_ms;
    CirculantRing<Zt, P1> m;
    Rp1 s[P1];
    Rp1_crt s_ms[P1];

    double noise = 0;
    double cumul_noise = 0;

    for (size_t idx = 0 ; idx < nb_tests ; ++idx)
    {
        m = CirculantRing<Zt, P1>::uniform_sample();
        for (size_t i = 0 ; i < P1 ; ++i)
        {
            s[i] = Rp1::sample_s(DENSITY_KEY);
            s_ms[i] = s[i];
        }
        rlwe.encrypt(s, m, Q, T, pow(2, var_in));
        //rlwe.encrypt(s, m, Q, T, 0);
        //rlwe.encrypt(s, m, Q, T, 32);
        //std::cout << "Noise before: " << lwe.noise(s, Q, T) << std::endl;

        rlwe_ms = rlwe.mod_switch<Rp1_crt, Qcrt, Q>();

        //std::cout << "Noise after: " << rlwe_ms.noise(s_ms, Qp, T) << std::endl;
        noise = rlwe_ms.noise(s_ms, Qcrt, T);
        cumul_noise += noise;
    }
    double var_out = log2(cumul_noise / nb_tests / P1);
    return var_out;
}

double exp_crt(const double var_in, const size_t nb_tests)
{
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

    for (size_t idx = 0 ; idx < nb_tests ; ++idx)
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

        c_p.encrypt(s_p, Xmp, Qcrt, T, pow(2, var_in));
        c_q.encrypt(s_q, Ymq, Qcrt, T, pow(2, var_in));
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
    double var_out = log2(cumul_noise / nb_tests / (P1*P2));
    return var_out;
}

double mod_switch_tensor(const double var_in, const size_t nb_tests)
{
    RLWE<Rp12_crt, 3> rlwe;
    Rp12LWE rlwe_ms;
    CirculantRing<Zt, P1*P2, FFT_DIM2> m;
    Rp12_crt s[3];
    Rp12 s_ms[3];

    double noise = 0;
    double cumul_noise = 0;

    for (size_t idx = 0 ; idx < nb_tests ; ++idx)
    {
        m = CirculantRing<Zt, P1*P2, FFT_DIM2>::uniform_sample();
        for (size_t i = 0 ; i < 3 ; ++i)
        {
            s[i] = Rp12_crt::sample_s(DENSITY_KEY);
            s_ms[i] = s[i];
        }
        rlwe.encrypt(s, m, Qcrt, T, pow(2, var_in));

        rlwe_ms = rlwe.mod_switch<Rp12, Qp, Qcrt>();

        noise = rlwe_ms.noise(s_ms, Qp, T);
        cumul_noise += noise;
    }
    double var_out = log2(cumul_noise / nb_tests / (P1*P2));
    return var_out;
}

double fun_exp_extract(const double var_in, const size_t nb_tests)
{
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

    for (size_t idx = 0 ; idx < nb_tests ; ++idx)
    {
        m = rand() % (P1*P2);
        for (size_t j = 0 ; j < P1*P2 ; ++j) Zm_coefs[j] = 0;
        Zm_coefs[m % (P1*P2)] = 1;
        CirculantRing<Zt, P1*P2, FFT_DIM2> Zm(Zm_coefs);

        rlwe_fun.encrypt(s_fun, Zm, Qp, T, pow(2, var_in));

        gen_funexpextract_key(S, s_fun, sp_fun, VARIANCE_FUNEXPEXTRACT);
        //std::cout << "Noise before: " << rlwe_fun.noise(s_fun, Qp, T) << std::endl;
        lwe_fun = fun_exp_extract(f, rlwe_fun, *S);

        //std::cout << "Noise after: " <<  lwe_fun.noise(sp_fun, Qp, T) << std::endl;
        noise = lwe_fun.noise(sp_fun, Qp, T);
        cumul_noise += noise;
    }
    delete S;
    delete[] Zm_coefs;
    free(f);
    double var_out = log2(cumul_noise / nb_tests);
    return var_out;
}

double sum_pow2(const double var_in, const size_t nb_tests)
{
    LWE c_i[6], res;
    CirculantRing<Zt, 1, 1> m_i[6];
    Rz s[P1];

    double noise = 0;
    double cumul_noise = 0;

    const size_t k = INPUT_BIT;
    int64_t coefs[k];
    if (k > 10)
    {
        for (size_t i = 0 ; i < k ; ++i)
            coefs[i] = 1;
    }
    else
    {
        coefs[0] = 1;
        for (size_t i = 1 ; i < k ; ++i)
            coefs[i] = 2 * coefs[i-1];
    }

    for (size_t idx = 0 ; idx < nb_tests ; ++idx)
    {
        for (size_t i = 0 ; i < P1 ; ++i)
            s[i] = Rz::sample_s(DENSITY_KEY);
        for (size_t i = 0 ; i < 6 ; ++i)
        {
            m_i[i] = rand() % 2;
            c_i[i].encrypt(s, m_i[i], Qp, T, pow(2, var_in));
        }

        res = combination(6, coefs, c_i);

        noise = res.noise(s, Qp, T);
        //noise = c_i[0].noise(s, Qp, T);
        cumul_noise += noise;
    }
    double var_out = log2(cumul_noise / nb_tests);
    return var_out;
}

double key_switch(const double var_in, const size_t nb_tests)
{
    RLWE<Rz, P1> in;
    double noise = 0;
    double cumul_noise = 0;

    Rz *s = new Rz[P1];
    for (size_t i = 0 ; i < P1 ; ++i)
        s[i] = Rz::sample_s(DENSITY_KEY);
    Rz *s2 = new Rz[N];
    for (size_t i = 0 ; i < N ; ++i)
        s2[i] = Rz::sample_s(DENSITY_KEY_SMALL);
    for (size_t idx = 0 ; idx < nb_tests ; ++idx)
    {
        KSWKeyLWE *S_lwe = new KSWKeyLWE(s, s2, 0);
        
        CirculantRing<Zt, 1, 1> m = CirculantRing<Zt, 1, 1>::uniform_sample();
        in.encrypt(s, m, Qp, T, pow(2, var_in));
        RLWE<Rz, N> out = in.key_switch(*S_lwe);

        noise = out.noise(s2, Qp, T);
        cumul_noise += noise;
        delete S_lwe;
    }
    delete[] s;
    delete[] s2;
    double var_out = log2(cumul_noise / nb_tests);
    return var_out;
}

double mod_switch_lwe(const double var_in, const size_t nb_tests)
{
    RLWE<Rz, N> Em_short;
    RLWE<CirculantRing<Zp12, 1, 1>, N> ab_short;
    CirculantRing<Zt, 1, 1> m;
    Rz s[N];
    CirculantRing<Zp12, 1, 1> s_ms[N];

    double noise = 0;
    double cumul_noise = 0;

    for (size_t idx = 0 ; idx < nb_tests ; ++idx)
    {
        m = CirculantRing<Zt, 1, 1>::uniform_sample();
        for (size_t i = 0 ; i < N ; ++i)
        {
            s[i] = Rz::sample_s(DENSITY_KEY);
            s_ms[i] = s[i];
        }
        Em_short.encrypt(s, m, Qp, T, pow(2, var_in));
        //std::cout << "Noise before: " << lwe.noise(s, Q, T) << std::endl;

        ab_short = Em_short.mod_switch<CirculantRing<Zp12, 1, 1>, P1*P2, Qp>();

        //std::cout << "Noise after: " << lwe_ms.noise(s_ms, P1*P2, T) << std::endl;
        noise = ab_short.noise(s_ms, P1*P2, T);
        cumul_noise += noise;
    }
    double var_out = log2(cumul_noise / nb_tests);
    return var_out;
}

typedef struct
{
    string name;
    double var_in;
    double var_out;
    size_t nb_tests;
    double (*function)(const double, const size_t);
} Step;


int main()
{
    Step steps[8] = {
        {    "ExtExpInner",  0.00, 38.50,  100, &ext_exp_inner},
        {   "ModSwitchAcc", 38.50,  6.41, 1000, &mod_switch_acc},
        {         "ExpCRT",  6.41, 24.82,  100, &exp_crt},
        {"ModSwitchTensor", 24.82, 60.82,  100, &mod_switch_tensor},
        {  "FunExpExtract", 60.82, 80.81,   50, &fun_exp_extract},
        {        "SumPow2", 80.81, 91.22, 1000, &sum_pow2},
        {      "KeySwitch", 91.22, 92.48,   50, &key_switch},
        {   "ModSwitchLWE", 92.48, 22.44, 1000, &mod_switch_lwe}
    };

    int seed = time(NULL);
    srand(seed);
    //std::cout << "seed " << seed << std::endl;
    print_param();
    Rp1GSW::G_init();

    double actual_var_out = 0.0;
    /*
    cout <<  "==== With expected variance as input" << endl;
    for (size_t i = 0 ; i < 8 ; ++i)
    {
        actual_var_out = steps[i].function(steps[i].var_in, steps[i].nb_tests);
        cout << steps[i].name << ":\tvar_in=2^" << steps[i].var_in << ",\tvar_out=2^" << actual_var_out << " (expected: 2^" << steps[i].var_out << "),\tnb_tests=" << steps[i].nb_tests << endl;
    }
    */
    cout <<  "==== With previous output variance as input" << endl;
    for (size_t i = 0 ; i < 8 ; ++i)
    {
        double previous_var_out = actual_var_out;
        actual_var_out = steps[i].function(actual_var_out, steps[i].nb_tests);
        cout << steps[i].name << ":\tvar_in=2^" << previous_var_out << ",\tvar_out=2^" << actual_var_out << " (expected: 2^" << steps[i].var_out << "),\tnb_tests=" << steps[i].nb_tests << endl;
    }

    return 0;
}
