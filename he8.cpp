#include <ctime>
#include <chrono>
#include <iostream>
#include <valgrind/callgrind.h>
#include "operations.h"
#include "ksw_key.h"
using namespace std;

int main(int argc, char *argv[])
{
    chrono::time_point<std::chrono::system_clock> start, end;
    chrono::duration<double> keygen_time, gate_time;
    //int pause;
    //cin >> pause;
    int seed = time(NULL);
    //seed = 1495639705;
    srand(seed);
    std::cerr << "seed " << seed << std::endl;
    cerr << "Generating gadgets... " << flush;
    Rp1GSW::G_init();
    Rp2GSW::G_init();
    cerr << "Done." << endl;
    //cin >> pause;

    const uint64_t q = Qp;
    const uint64_t t = T;
    const size_t k = INPUT_BIT;
    const size_t nb_iter = (argc > 1) ? atoi(argv[1]) : 5;

    /** Input data **/
    LWE c_i[k];
    Fun F;
    //std::cout << "F(" << m << ") = F(" << m*(P1*P2/t) << ") = " << F(m*(P1*P2/t)) << std::endl;
    Rp12 Rf(F);
    fftw_complex *f = align_alloc<fftw_complex>(ALIGNMENT, FFT_DIM2);
    Rf.compute_fft(f);

    int64_t coefs[k];
    if (k > 10)
    {
        cerr << "Computing 63 bits threshold gate" << endl;
        for (size_t i = 0 ; i < k ; ++i)
            coefs[i] = 1;
    }
    else
    {
        cerr << "Computing 6 bits parity gate" << endl;
        coefs[0] = 1;
        for (size_t i = 1 ; i < k ; ++i)
            coefs[i] = 2 * coefs[i-1];
    }

    start = chrono::system_clock::now();
    /** Generate secret keys **/
    cerr << "Generating secret keys... " << flush;
    Rz *s = align_alloc<Rz>(ALIGNMENT, P1);
    for (size_t i = 0 ; i < P1 ; ++i)
        s[i] = Rz::sample_s(DENSITY_KEY);
    Rz *s2 = align_alloc<Rz>(ALIGNMENT, N);
    for (size_t i = 0 ; i < N ; ++i)
        s2[i] = Rz::sample_s(DENSITY_KEY_SMALL);
    Rp1 s_p = Rp1::sample_s(DENSITY_KEY_ACC);
    Rp2 s_q = Rp2::sample_s(DENSITY_KEY_ACC);
    Rp12 s_pq[3];
    crt_key(s_pq, s_p, s_q);
    cerr << "Done." << endl;
    //cout << "s=("; for (size_t i = 0 ; i < P1 ; ++i) cout << s[i].get_data()[0] << ","; cout << ")" << endl;
    //cin >> pause;

    /** Derive key material **/
    cerr << "Generating LWE key-switching key... " << flush;
    KSWKeyLWE *S_lwe = new KSWKeyLWE(s, s2, VARIANCE_KSWLWE);
    cerr << "Done." << endl;
    //cin >> pause;

    cerr << "Generating bootstrapping keys... " << flush;
    Rp1GSW *Xsi = new Rp1GSW[N];
    gen_bootstrapping_keys<Rp1, N, P1>(q, t, Xsi, s2, s_p, VARIANCE_ACC);
    cerr << "1/2 " << flush;
    //cin >> pause;
    Rp2GSW *Ysi = new Rp2GSW[N];
    gen_bootstrapping_keys<Rp2, N, P2>(q, t, Ysi, s2, s_q, VARIANCE_ACC);
    cerr << "Done." << endl;
    free(s2);
    //cin >> pause;

    cerr << "Generating key-switching keys... " << flush;
    KSWKeyRp1 *KSp1 = new KSWKeyRp1[P1];
    gen_keyswitching_keys<Rp1, P1>(KSp1, s_p, VARIANCE_ACC);
    cerr << "1/2 " << flush;
    //cin >> pause;
    KSWKeyRp2 *KSp2 = new KSWKeyRp2[P2];
    gen_keyswitching_keys<Rp2, P2>(KSp2, s_q, VARIANCE_ACC);
    cerr << "Done." << endl;
    //cin >> pause;

    cerr << "Generating function extraction key... " << flush;
    KSWKeyRp12 *S = new KSWKeyRp12();
    gen_funexpextract_key(S, s_pq, s, VARIANCE_FUNEXPEXTRACT);
    cerr << "Done." << endl;
    //cin >> pause;
    end = chrono::system_clock::now();
    keygen_time = end - start;
    cerr << "keygen time: " << keygen_time.count() << "s" << endl;


    size_t nb_success = 0;
    chrono::duration<double> total_time = chrono::duration<double>::zero();
    for (size_t iter = 1 ; iter <= nb_iter ; ++iter)
    {
        cerr << "   =======    Test " << iter << "/" << nb_iter << "    ========" << endl;
        /** Setup values **/
        cerr << "Generating values... " << flush;
        CirculantRing<Zt, 1, 1> m_i[k];
        for (size_t i = 0 ; i < k ; ++i)
        {
            m_i[i] = rand() % 2;
            c_i[i].encrypt(s, m_i[i], q, t, VARIANCE_INPUT);
            //std::cout << m_i[i] << " "; //<< " (" << c_i[i] << ") ";
        }
        int64_t m = 0;
        for (size_t i = 0 ; i < k ; ++i)
            m += coefs[i] * m_i[i].get_data()[0].get_value();
        cerr << "Done." << endl;
        cerr << "Computing F(" << m << ") -> " << F.eval(m) << endl;
        cerr << "Timings are:" << endl;

        /** Launch gate **/
        CALLGRIND_START_INSTRUMENTATION;
        CALLGRIND_TOGGLE_COLLECT;
        start = chrono::system_clock::now();
        LWE result = gate(k, coefs, c_i, f, *S_lwe, Xsi, KSp1, Ysi, KSp2, *S);
        end = chrono::system_clock::now();
        gate_time = end - start;
        CALLGRIND_TOGGLE_COLLECT;
        CALLGRIND_STOP_INSTRUMENTATION;

        CirculantRing<Zt, 1, 1> d;
        result.decrypt(s, d, q, t);

        total_time += gate_time;
        cerr << "Total gate time " << gate_time.count() << "s" << endl;
        cout << "Actual result F(" << m << ") = " << d.get_data()[0] << " (expected " << F.eval(m) << ")" << endl << endl;
        if ((d - CirculantRing<Zt, 1, 1>(F.eval(m))).is_zero())
            ++nb_success;
    }

    cout << endl << "Success: " << nb_success << "/" << nb_iter << " (" << nb_success * 100 / nb_iter << " %), average gate time: " << total_time.count()/nb_iter << "s" <<  endl << endl;


    free(f);
    free(s);
    delete S;
    delete[] Xsi;
    delete[] Ysi;
    delete[] KSp1;
    delete[] KSp2;
    delete S_lwe;

    return 0;
}

