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

    double EPS = 5;

    // Хард код, параметры s2 на первом наблюдении
    three_stars[0] = -1.34542432318288867188e+13 - -1.34542432318288867188e+13/ EPS;
    three_stars[1] = 2.74369376973100244141e+12 + 2.74369376973100244141e+12 / EPS;
    three_stars[2] = 1.17902665512916289062e+13 - 1.17902665512916289062e+13 / EPS;
    three_stars[3] = 9.63030905550202805898e+03 + 9.63030905550202805898e+03 / EPS;
    three_stars[4] = 2.38743901750857585284e+04 - 2.38743901750857585284e+04 / EPS;

    // Хард код, параметры s38 на первом наблюдении
    three_stars[5] = 4.0707e+12 - 4.0707e+12/ EPS;
    three_stars[6] = 3.11869e+13 + 3.11869e+13 / EPS;
    three_stars[7] = 2.54138e+12 - 2.54138e+12 / EPS;
    three_stars[8] = 18926.6 + 18926.6 / EPS;
    three_stars[9] = -2617.61 + -2617.61 / EPS;

    // Хард код, параметры s55 на первом наблюдении
    three_stars[10] = 3.06247e+13 - 3.06247e+13/ EPS;
    three_stars[11] = -3.22222e+12 -3.22222e+12 / EPS;
    three_stars[12] = 1.69223e+13 - 1.69223e+13 / EPS;
    three_stars[13] = 1642.9 + 1642.9 / EPS;
    three_stars[14] =-16531.6  -16531.6  / EPS;
    three_stars[15] = 1e32;

    gauss_newton(three_stars, 3);
}