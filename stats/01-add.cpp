#include "../fft.h"
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

    LWE c_i[6], res;
    CirculantRing<Zt, 1, 1> m_i[6];
    Rz s[P1];

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

    double noise = 0;
    double cumul_noise = 0;

    for (size_t idx = 0 ; idx < NB_TESTS ; ++idx)
    {
        for (size_t i = 0 ; i < P1 ; ++i)
            s[i] = Rz::sample_s(DENSITY_KEY);
        for (size_t i = 0 ; i < 6 ; ++i)
        {
            m_i[i] = rand() % 2;
            c_i[i].encrypt(s, m_i[i], Qp, T, std::pow(2, 80.81));
            //c_i[i].encrypt(s, m_i[i], Qp, T, 1024);
            //c_i[i].encrypt(s, m_i[i], Qp, T, 256.0);
            //std::cout << m_i[i] << " "; //<< " (" << c_i[i] << ") ";
        }
        //std::cout << "Noise before: " << lwe1.noise(s, Q, T) << std::endl;
        //std::cout << "Noise before: " << lwe2.noise(s, Q, T) << std::endl;

        res = combination(6, coefs, c_i);

        noise = res.noise(s, Qp, T);
        //noise = c_i[0].noise(s, Qp, T);
        cumul_noise += noise;
        //std::cout << noise << std::endl;
        //std::cout << cumul_noise << std::endl;
    }
    std::cout << cumul_noise << " " << cumul_noise / NB_TESTS << " = 2^" << std::log2(cumul_noise / NB_TESTS) << std::endl;

    return 0;
}
