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
    double three_stars[19];

    double EPS = 5;

    // Хард код, параметры s2 на первом наблюдении
    three_stars[0] = -1.34542432318288867188e+13 - -1.34542432318288867188e+13/ EPS;
    three_stars[1] = 2.74369376973100244141e+12 + 2.74369376973100244141e+12 / EPS;
    three_stars[2] = 1.17902665512916289062e+13 - 1.17902665512916289062e+13 / EPS;
    three_stars[3] = 9.63030905550202805898e+03 + 9.63030905550202805898e+03 / EPS;
    three_stars[4] = 2.38743901750857585284e+04 - 2.38743901750857585284e+04 / EPS;
    three_stars[5] = 5.660631180e+03 - 5.660631180e+03 / EPS;

    // Хард код, параметры s38 на первом наблюдении
    three_stars[6] = 4.0707e+12 - 4.0707e+12/ EPS;
    three_stars[7] = 3.11869e+13 + 3.11869e+13 / EPS;
    three_stars[8] = 2.54138e+12 - 2.54138e+12 / EPS;
    three_stars[9] = 18926.6 + 18926.6 / EPS;
    three_stars[10] = -2617.61 + -2617.61 / EPS;
    three_stars[11] = 4412.44 - 4412.44/ EPS;

    // Хард код, параметры s55 на первом наблюдении
    three_stars[12] = 3.06247e+13 - 3.06247e+13/ EPS;
    three_stars[13] = -3.22222e+12 -3.22222e+12 / EPS;
    three_stars[14] = 1.69223e+13 - 1.69223e+13 / EPS;
    three_stars[15] = 1642.9 + 1642.9 / EPS;
    three_stars[16] =-16531.6  -16531.6  / EPS;
    three_stars[17] = -7379.3  + 7379.3 / EPS;
    three_stars[18] = 1e32;

    int x, y, z, vx, vy, vz;

    scanf("%d %d %d %d %d %d", &x, &y, &z, &vx, &vy, &vz);

    int conditions[6] = {x, y, z, vx, vy, vz};
    int size=0;
    for(int i=0; i<6; i++)
        if(conditions[i]==1) size++;

    int num_star = 1;

    double result_vector[size*num_star+1];

    int tmp = 0;

    for(int i=0; i<6*num_star; i++)
    {
        if(conditions[i % 6] == 1)
            result_vector[tmp++] = three_stars[i];
    }
    result_vector[size*num_star] = 1e32;

    gauss_newton(result_vector, num_star, conditions);
}