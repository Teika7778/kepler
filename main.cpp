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
double three_stars[16];

    // Хард код, параметры s2 на первом наблюдении
    three_stars[0] = 2.14992e+13;
    three_stars[1] = 3.59558e+13;
    three_stars[2] = 3.17322e+12;
    three_stars[3] = 3.96405e+06;
    three_stars[4] = 1.91863e+06;

    // Хард код, параметры s38 на первом наблюдении
    three_stars[5] = 8.03064e+13 ;
    three_stars[6] = -8.03314e+13;
    three_stars[7] = 1.52503e+13 ;
    three_stars[8] = 365942;
    three_stars[9] = -2.54613e+06;

    // Хард код, параметры s55 на первом наблюдении
    three_stars[10] = -1.94304e+14;
    three_stars[11] = 7.44603e+13;
    three_stars[12] = -8.00027e+13;
    three_stars[13] = 478498;
    three_stars[14] = 567817;
    three_stars[15] = 8e36;

    gauss_newton(three_stars, 3);
}