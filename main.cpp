#include <cmath>
#include <iostream>

#include "constans.hpp"
#include "struct.hpp"
#include "normalize.hpp"
#include "transform.hpp"
#include "diff.hpp"
#include "gauss_newton.hpp"
#include "integration.hpp"

int main()
{
    double parameters[7];

    // Хард код, параметры s2 на первом наблюдении
    parameters[0] = -1.34542432318288867188e+13;
    parameters[1] = 2.74369376973100244141e+12;
    parameters[2] = 1.17902665512916289062e+13;
    parameters[3] = 9.63030905550202805898e+03;
    parameters[4] = 2.38743901750857585284e+04;
    parameters[5] = 5.6606311809e+03;
    parameters[6] = 1e32;

    gauss_newton(parameters);
}