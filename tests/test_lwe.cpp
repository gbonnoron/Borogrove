#include <iostream>
#include <ctime>
#include <cmath>
#include "../ksw_key.h"
#include "../gadget.h"
#include "../circulant_ring.h"
#include "../predicate.h"
#include "../operations.h"

int main()
{
    int seed = time(NULL);
    srand(seed);
    //std::cout << "seed " << seed << std::endl;
    Rz s[P1];

    for (size_t i = 0 ; i < P1 ; ++i)
        s[i] = Rz::sample_s();

    //std::cout << "s "; for (size_t i = 0 ; i < P1 ; ++i) std::cout << s[i] << " "; std::cout << std::endl;

    for (unsigned int i = 0 ; i < 3 ; ++i)
    {
        //std::cout << "Test " << i << std::endl;

        /** Testing basic encryption/decryption **/
        CirculantRing<Zt, 1, 1> m = CirculantRing<Zt, 1, 1>::uniform_sample();
        //std::cout << "m\t" << m << std::endl;
        LWE lwe;
        lwe.encrypt(s, m, Qp, T, VARIANCE_INPUT);
        CirculantRing<Zt, 1, 1> d;
        lwe.decrypt(s, d, Qp, T);
        //std::cout << "d\t" << d << std::endl;
        if (! (m - d).is_zero())
        {
            std::cerr << "Encrypt/decrypt failed." << std::endl;
            return 1;
        }

        /** Testing addition **/
        CirculantRing<Zt, 1, 1> m1 = CirculantRing<Zt, 1, 1>::uniform_sample();
        //std::cout << "m1\t" << m1 << std::endl;
        LWE lwe1;
        lwe1.encrypt(s, m1, Qp, T, VARIANCE_INPUT);
        lwe.decrypt(s, d, Qp, T);
        //std::cout << "d1\t" << d << std::endl;

        lwe1 += lwe;
        lwe1.decrypt(s, d, Qp, T);
        //std::cout << "m+m1\t" << d << std::endl;
        if (!(d - (m + m1)%T).is_zero())
        {
            std::cerr << "Addition failed." << std::endl;
            return 1;
        }

        /** Testing modulus switching **/
        Rz s_ms2[N];
        CirculantRing<Zp12, 1, 1> s_ms[N];
        for (size_t i = 0 ; i < N ; ++i)
        {
            s_ms2[i] = Rz::sample_s();
            s_ms[i] = s_ms2[i];
            //std::cout << s[i] << " -> " << s_ms[i] << std::endl;
        }
        RLWE<Rz, N> lwe2;
        RLWE<CirculantRing<Zp12, 1, 1>, N> lwe_ms;
        lwe2.encrypt(s_ms2, m, Qp, T, VARIANCE_INPUT);
        lwe_ms = lwe2.mod_switch<CirculantRing<Zp12, 1, 1>, P1*P2, Qp>();
        //std::cout << lwe << std::endl << std::endl;
        //std::cout << lwe2 << std::endl;
        lwe_ms.decrypt(s_ms, d, P1*P2, T);
        //std::cout << "mm\t" << d << std::endl;
        //std::cout << "noise before " << std::log2(lwe.noise(s, Qp, T)) << std::endl;
        //std::cout << "noise after " << std::log2(lwe2.noise(s_ms, P1*P2, T)) << std::endl;
        if (! (m - d).is_zero())
        {
            std::cerr << "ModSwitch failed." << std::endl;
            return 1;
        }

        /** Testing key switching **/
        Rp12LWE lwe3;
        RLWE<Rp12, 1> lwe4;
        Rp12 s1[3];
        for (size_t j=0 ; j<3 ; ++j)
        {
            s1[j] = Rp12::sample_s();
            //s1[j].decomp_fft();
        }
        Rp12 s2[1];
        s2[0] = Rp12::sample_s();
        //std::cout << "s1_0 " << s1[0] << std::endl;
        //std::cout << "s1_1 " << s1[1] << std::endl;
        //std::cout << "s1_2 " << s1[2] << std::endl;
        //std::cout << "s2_0 " << s2[0] << std::endl;
        CirculantRing<Zt, P1*P2, FFT_DIM2>  m_ksw = CirculantRing<Zt, P1*P2, FFT_DIM2>::uniform_sample();
        lwe3.encrypt(s1, m_ksw, Qp, T, VARIANCE_INPUT);

        KSWKeyRp12 *ksw = new KSWKeyRp12(s1, s2, VARIANCE_FUNEXPEXTRACT);
        lwe4 = lwe3.key_switch(*ksw);
        delete ksw;
        CirculantRing<Zt, P1*P2, FFT_DIM2> d_ksw;
        lwe4.decrypt(s2, d_ksw, Qp, T);
        //std::cout << "m\t" << m_ksw << std::endl;
        //std::cout << "d\t" << d_ksw << std::endl;
        //std::cout << "noise " << log2(lwe4.noise(s2, Qp, T)) << std::endl;
        if (! (m_ksw - d_ksw).is_zero())
        {
            std::cerr << "KeySwitch failed." << std::endl;
            return 1;
        }

        /** Testing Galois **/
        uint64_t alpha = 2;
        //std::cout << "alpha " << alpha << std::endl;
        CirculantRing<Zt, P1> m_galois = CirculantRing<Zt, P1>::uniform_sample();
        Rp1 s_g[1];
        s_g[0] = Rp1::sample_s();
        Rp1 s_galois[1] = {s_g[0].galois(alpha)};
        Rp1LWE lwe_galois;
        lwe_galois.encrypt(s_g, m_galois, Q, T, VARIANCE_INPUT);
        lwe_galois.galois_inplace(alpha);
        m_galois.galois_inplace(alpha);
        CirculantRing<Zt, P1> d_galois;
        lwe_galois.decrypt(s_galois, d_galois, Q, T);
        if (! (d_galois - m_galois).is_zero())
        {
            std::cerr << "Galois failed." << std::endl;
            return 1;
        }
        /**/

        size_t m4 = rand() % T;
        size_t m3 = (double)m4 * (double)(P1*P2) / (double)T;
        //std::cout << "m3 " << m3 << std::endl;
        Zt *Zm_coefs = new Zt[P1*P2];
        for (size_t j = 0 ; j < P1*P2 ; ++j) Zm_coefs[j] = 0;
        Zm_coefs[m3 % (P1*P2)] = 1;
        CirculantRing<Zt, P1*P2, FFT_DIM2> Zm(Zm_coefs, false);
        //std::cout << "Z^m\t" << Zm << std::endl;

        /** Testing ExpCRT **/
        Zt Xmp_coefs[P1];
        for (size_t j = 0 ; j < P1 ; ++j) Xmp_coefs[j] = 0;
        Xmp_coefs[m3 % P1] = 1;
        CirculantRing<Zt, P1> Xmp(Xmp_coefs);
        //std::cout << "X^mp\t" << Xmp << std::endl;
        Zt Ymq_coefs[P2];
        for (size_t j = 0 ; j < P2 ; ++j) Ymq_coefs[j] = 0;
        Ymq_coefs[m3 % P2] = 1;
        CirculantRing<Zt, P2> Ymq(Ymq_coefs);
        //std::cout << "Y^mq\t" << Ymq << std::endl;

        Rp1 s_p = Rp1::sample_s();
        //std::cout << "s_p " << s_p << std::endl;
        Rp2 s_q = Rp2::sample_s();
        //std::cout << "s_q " << s_q << std::endl;
        Rp12 s_pq[3];
        crt_key(s_pq, s_p, s_q);
        Rp1LWE c_p;
        c_p.encrypt(&s_p, Xmp, Q, T, VARIANCE_INPUT);
        RLWE<Rp1_crt, 1> c_p_crt = c_p.template mod_switch<Rp1_crt, Qcrt, Q>();
        Rp2LWE c_q;
        c_q.encrypt(&s_q, Ymq, Q, T, VARIANCE_INPUT);
        RLWE<Rp2_crt, 1> c_q_crt = c_q.template mod_switch<Rp2_crt, Qcrt, Q>();
        RLWE<Rp12_crt, 3> result;
        exp_crt(result, c_p_crt, c_q_crt);
        Rp12LWE c_pq = result.template mod_switch<Rp12, Qp, Qcrt>();
        //std::cout << "s_pq " << s_pq[0] << ", " << s_pq[1] << ", " << s_pq[2] << std::endl;
        CirculantRing<Zt, P1*P2, FFT_DIM2> dec;
        c_pq.decrypt(s_pq, dec, Qp, T);
        //std::cout << "crt noise " << std::log2(c_pq.noise(s_pq, Qp, T)) << std::endl;
        //std::cout << "Z^m\t" << dec << std::endl;
        if (! (Zm - dec).is_zero())
        {
            std::cerr << "ExpCRT failed." << std::endl;
            return 1;
        }

        /** Testing function extraction **/
        Fun F;
        Rp12 Rf(F);
        fftw_complex *f = align_alloc<fftw_complex>(ALIGNMENT, FFT_DIM2);
        Rf.compute_fft(f);
        Rp12LWE rlwe_fun;
        Rp12 s_fun[3];
        for (size_t j=0 ; j<3 ; ++j)
        {
            s_fun[j] = Rp12::sample_s();
            //std::cout << s_fun[j] << std::endl;
        }
        rlwe_fun.encrypt(s_fun, Zm, Qp, T, VARIANCE_INPUT);
        //rlwe_fun.decrypt(s_fun, dec, Qp, T);
        //std::cout << "dec " << dec << std::endl;

        Rz sp_fun[P1];
        for (size_t j = 0 ; j < P1 ; ++j)
        {
            sp_fun[j] = Rz::sample_s();
            //std::cout << sp_fun[j] << " ";
        }
        //std::cout << std::endl;
        Rp12 spp_fun[1];
        Zqp *coefs = new Zqp[P1*P2];
        for (size_t k = 0 ; k < P1*P2 ; ++k)
            coefs[k] = 0;
        for (size_t k = 0 ; k < P1 ; ++k)
            coefs[k * P2] = sp_fun[k].get_data()[0];
        spp_fun[0] = Rp12(coefs, false);
        KSWKeyRp12 *S = new KSWKeyRp12();
        //gen_funexpextract_key(S, s_fun, sp_fun, FUNEXPEXTRACT);
        //LWE lwe_fun = fun_exp_extract(f, rlwe_fun, *S);
        //std::cout << "noise before\t" << std::log2(rlwe_fun.noise(s_fun, Qp, T)) << std::endl;
        gen_funexpextract_key(S, s_pq, sp_fun, VARIANCE_FUNEXPEXTRACT);
        LWE lwe_fun = fun_exp_extract(f, c_pq, *S);
        free(f);
        delete S;

        CirculantRing<Zt, 1, 1> dec_fun;
        lwe_fun.decrypt(sp_fun, dec_fun, Qp, T);

        //std::cout << "m\t" << m4 << std::endl;
        //std::cout << "d_exp\t" << F.eval(m4) << std::endl;
        //std::cout << "d_fun\t" << dec_fun << std::endl;
        //std::cout << "noise after\t" << std::log2(lwe_fun.noise(sp_fun, Qp, T)) << std::endl;
        if (! (dec_fun - CirculantRing<Zt, 1, 1>(F.eval(m4))).is_zero())
        {
            std::cerr << "FunExpExtract failed." << std::endl;
            return 1;
        }
        /**/
    }

    return 0;
}
