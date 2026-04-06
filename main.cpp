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

    double EPS = 5000000000;

    // Хард код, параметры s2 на первом наблюдении
    three_stars[0] = (2.14992e+13 - 2.14992e+13/ EPS) / AU;
    three_stars[1] = (3.59558e+13 + 3.59558e+13 / EPS) / AU;
    three_stars[2] = (3.17322e+12 - 3.17322e+12 / EPS) / AU;
    three_stars[3] = (3.96405e+06 + 3.96405e+06 / EPS) / AU * DAY * YEAR;
    three_stars[4] = (1.91863e+06 - 1.91863e+06 / EPS) / AU * DAY * YEAR;
    three_stars[5] = (-1.98566e+06) / AU * DAY * YEAR;

    // Хард код, параметры s38 на первом наблюдении
    three_stars[6] = (8.03064e+13 - 8.03064e+13/ EPS ) / AU;
    three_stars[7] = (-8.03314e+13 -8.03314e+13 / EPS) / AU;
    three_stars[8] = (1.52503e+13 - 1.52503e+13 / EPS) / AU;
    three_stars[9] = (365942 + 365942 / EPS) / AU * DAY * YEAR;
    three_stars[10] = (-2.54613e+06 + -2.54613e+06 / EPS) / AU * DAY * YEAR;
    three_stars[11] = (-39845.1) / AU * DAY * YEAR;

    // Хард код, параметры s55 на первом наблюдении
    three_stars[12] = (-1.94304e+14 + -1.94304e+14/ EPS) / AU;
    three_stars[13] = (7.44603e+13 -7.44603e+13 / EPS) / AU;
    three_stars[14] = (-8.00027e+13 + -8.00027e+13 / EPS) / AU;
    three_stars[15] = (478498 + 478498 / EPS) / AU * DAY * YEAR;
    three_stars[16] = (567817  -567817  / EPS) / AU * DAY * YEAR;
    three_stars[17] = (577556) / AU * DAY * YEAR;
    three_stars[18] = 8e36 / M_SUN;

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
    result_vector[size*num_star] = 8e36 / M_SUN;

    gauss_newton(result_vector, num_star, conditions);
}