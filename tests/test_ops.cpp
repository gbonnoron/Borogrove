#include <iostream>
#include <ctime>
#include <cstring>

#include "../operations.h"

typedef CirculantRing<Zt, P1> MyRpt;

int main()
{
    int seed = time(NULL);
    srand(seed);
    //std::cout << "seed " << seed << std::endl;
    Rp1GSW::G_init();

    Rp1 s[1];
    s[0] = Rp1::sample_s();
    //std::cout << "s " << s[0] << std::endl;

    for (unsigned int idx = 0 ; idx < 5 ; ++idx)
    {
        //std::cout << "Test " << idx << std::endl;

        /** Testing basic encryption and decryption **/
        size_t m = rand() % (P1);
        Zt Xm_coefs[P1];
        for (size_t j = 0 ; j < P1 ; ++j) Xm_coefs[j] = 0;
        Xm_coefs[m] = 1;
        MyRpt Xm(Xm_coefs);
        //std::cout << "Xm\t" << Xm << std::endl;
        Rp1LWE rlwe;
        rlwe.encrypt(s, Xm, Q, T, VARIANCE_INPUT);
        MyRpt d;
        rlwe.decrypt(s, d, Q, T);
        //std::cout << "Xm\t" << d << std::endl;
        if (! (Xm - d).is_zero())
        {
            std::cerr << "RLWE Encrypt/decrypt failed." << std::endl;
            return 1;
        }

        size_t mp = rand() % (P1);
        Zt Xmp_coefs[P1];
        for (size_t j = 0 ; j < P1 ; ++j) Xmp_coefs[j] = 0;
        Xmp_coefs[mp] = 1;
        MyRpt Xmp(Xmp_coefs);
        //std::cout << "Xmp\t" << Xmp << std::endl;
        Rp1GSW rgsw;
        rgsw.encrypt(s, Xmp, VARIANCE_INPUT);
        rgsw.decrypt(s, d);
        //std::cout << "Xmp\t" << d << std::endl;
        if (! (Xmp - d).is_zero())
        {
            std::cerr << "RGSW Encrypt/decrypt failed." << std::endl;
            return 1;
        }

        /** Testing ExtMult **/
        RLWE<Rp1, 1> rlwe2 = rlwe.ext_mult(rgsw);
        Zt Xm_mp_coefs[P1];
        for (size_t j = 0 ; j < P1 ; ++j) Xm_mp_coefs[j] = 0;
        Xm_mp_coefs[(m + mp) % P1] = 1;
        MyRpt Xm_mp(Xm_mp_coefs);
        rlwe2.decrypt(s, d, Q, T);
        //std::cout << "noise extmult " << std::log2(rlwe2.noise(s, Q, T)) << std::endl;
        //std::cout << "Xm+mp\t" << Xm_mp << std::endl;
        //std::cout << "Xm+mp\t" << d << std::endl;
        if (! (Xm_mp - d).is_zero())
        {
            std::cerr << "ExtMult failed." << std::endl;
            return 1;
        }

        /** Testing ExtExpMultAdd **/
        uint64_t alpha = 2;
        uint64_t beta = (size_t)Zp1(alpha).inv();
        //std::cout << alpha << " " << beta << std::endl;
        Rp1 s_alpha[1], s_beta[1];
        s_alpha[0] = s[0].galois(alpha);
        s_beta[0]  = s[0].galois(beta);
        KSWKeyRp1 S_alpha(s_alpha, s, VARIANCE_ACC);
        KSWKeyRp1 S_beta(s_beta, s, VARIANCE_ACC);

        rlwe.ext_exp_mult_add(rgsw, alpha, S_alpha, beta, S_beta);
        Zt Xam_mp_coefs[P1];
        for (size_t j = 0 ; j < P1 ; ++j) Xam_mp_coefs[j] = 0;
        Xam_mp_coefs[(alpha * mp + m) % P1] = 1;
        MyRpt Xam_mp(Xam_mp_coefs);
        rlwe.decrypt(s, d, Q, T);
        //std::cout << "Xam+mp\t" << Xam_mp << std::endl;
        //std::cout << "Xam+mp\t" << d << std::endl;
        if (! (Xam_mp - d).is_zero())
        {
            std::cerr << "ExtExpMultAdd failed." << std::endl;
            return 1;
        }

        /** Testing ExtExpInner **/
        KSWKeyRp1 S[P1];
        gen_keyswitching_keys<Rp1, P1>(S, s[0], VARIANCE_ACC);

        Rp1GSW C[N];
        Zqp x[N];
        Rz s_x[N];
        Zp1 x2[N];
        for (size_t i = 0 ; i < N ; ++i)
        {
            x[i] = rand();
            x2[i] = x[i] % P1;
            s_x[i] = x[i];
        }
        gen_bootstrapping_keys<Rp1, N, P1>(Q, T, C, s_x, s[0], VARIANCE_ACC);
        Zp1 y[N];
        for (size_t i = 0 ; i < N ; ++i)
            y[i] = rand();

        //std::cout << "x "; for (size_t i = 0 ; i < N ; ++i) std::cout << x2[i] << " "; std::cout << std::endl;
        //std::cout << "y "; for (size_t i = 0 ; i < N ; ++i) std::cout << y[i] << " "; std::cout << std::endl;
        Zp1 xy = 0;
        for (size_t i = 0 ; i < N ; ++i)
            xy += x2[i] * y[i];
        //std::cout << "xy " << xy << std::endl;
        Zt sol_coefs[P1];
        for (size_t j = 0 ; j < P1 ; ++j) sol_coefs[j] = 0;
        sol_coefs[(size_t)xy] = 1;
        MyRpt sol(sol_coefs);
        //std::cout << "Xxy\t" << sol << std::endl;

        Rp1LWE c = ext_exp_inner<Rp1, Zp1, P1>(Q, T, N, y, C, S);
        c.decrypt(s, d, Q, T);
        //std::cout << "c\t" << d << std::endl;
        if (! (sol - d).is_zero())
        {
            std::cerr << "ExtExpInner failed." << std::endl;
            return 1;
        }
        /**/
    }

    return 0;
}
