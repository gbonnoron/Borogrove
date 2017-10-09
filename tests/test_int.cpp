#include <iostream>
#include "../integer_mod.h"

int main ()
{
    Zq a(10);
    Zq b(3);
    Zq c(0);
    Zq d(32);
    Zq e(-32);
    Zq eh, el;
    int q2 = 11;

    std::cout << "a " << a << std::endl;
    std::cout << "b " << b << std::endl;
    std::cout << "c " << c << std::endl;
    std::cout << "d " << d << std::endl;
    std::cout << "e " << e << std::endl;
    std::cout << "a-b " << a-b << std::endl;
    std::cout << "b-a " << b-a << std::endl;
    std::cout << "a * b " << a*b << std::endl;
    std::cout << "a % q2 " << a%q2 << std::endl;
    std::cout << "b % q2 " << b%q2 << std::endl;
    Zq res = a * 2;
    Zq res0 = b * 2;
    std::cout << "a * 2 " << a*2 << std::endl;
    std::cout << "b * 2 " << b*2 << std::endl;
    std::cout << "res " << res << std::endl;
    std::cout << "res0 " << res0 << std::endl;
    Zq res1 = a / 2;
    Zq res2 = b / 2;
    std::cout << "a / 2 " << a/2 << std::endl;
    std::cout << "b / 2 " << b/2 << std::endl;
    std::cout << "res1 " << res1 << std::endl;
    std::cout << "res2 " << res2 << std::endl;
    Zq res3 = (a/2) % 5;
    Zq res4 = (b/2) % 5;
    std::cout << "(a / 2) % 5 " << (a/2)%5 << std::endl;
    std::cout << "(b / 2) % 5 " << (b/2)%5 << std::endl;
    std::cout << "res3 " << res3 << std::endl;
    std::cout << "res4 " << res4 << std::endl;

    int64_t B = e.split(eh, el);
    std::cout << e << " = " << eh << " * " << B << " + " << el << std::endl;

    return 0;
}
