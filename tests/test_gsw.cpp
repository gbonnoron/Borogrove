#include <iostream>
#include <ctime>

#include "../circulant_ring.h"
#include "../param.h"
#include "../gadget.h"
#include "../rgsw.h"

typedef CirculantRing<Zt, P1> MyRpt;

int main()
{
    int seed = time(NULL);
    srand(seed);
    //std::cout << "seed " << seed << std::endl;
    Rp1GSW::G_init();

    Rp1 s[N];
    for (size_t i = 0 ; i < N ; ++i)
        s[i] = Rp1::sample_s();

    //std::cout << "s "; for (size_t i = 0 ; i < N ; ++i) std::cout << s[i]; std::cout << std::endl;

    for (unsigned int i = 0 ; i < 10 ; ++i)
    {
        //std::cout << "Test " << i << std::endl;

        /** Testing basic encryption/decryption **/
        MyRpt m = MyRpt::uniform_sample();
        //std::cout << "m\t" << m << std::endl;
        Rp1GSW rgsw;
        rgsw.encrypt(s, m, VARIANCE_INPUT);
        MyRpt d;
        rgsw.decrypt(s, d);
        //std::cout << "d\t" << d << std::endl;
        if (! (m - d).is_zero())
        {
            std::cerr << "Encrypt/decrypt failed." << std::endl;
            return 1;
        }
    }

    return 0;
}
