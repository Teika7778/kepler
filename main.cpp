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

    // Хард код, параметры s2 на первом наблюдении
    three_stars[0] = 0.126;  // a
    three_stars[1] = 0.884;  // e
    three_stars[2] = 71.36;   // w
    three_stars[3] = 234.50;  // omega
    three_stars[4] = 136.78;  // i
    three_stars[5] = 2002.32; // T0

    // Хард код, параметры s38 на первом наблюдении
    three_stars[6] = 0.140;  // a
    three_stars[7] = 0.818;  // e
    three_stars[8] = 18.4;   // w
    three_stars[9] = 101.8;  // omega
    three_stars[10] = 166.22;  // i
    three_stars[11] = 2003.30; // T0

    // Хард код, параметры s55 на первом наблюдении
    three_stars[12] = 0.109;  // a
    three_stars[13] = 0.74;  // e
    three_stars[14] = 133.5;   // w
    three_stars[15] = 129.9;  // omega
    three_stars[16] = 141.7;  // i
    three_stars[17] = 2009.31; // T0
    three_stars[18] = 8e36;

    int a, e, w, omega, i, T0;

    scanf("%d %d %d %d %d %d", &a, &e, &w, &omega, &i, &T0);

    int conditions[6] = {a, e, w, omega, i, T0};
    int size=0;
    for(int i=0; i<6; i++)
        if(conditions[i]==1) size++;

    int num_star = 3;

    double result_vector[size*num_star+1];

    int tmp = 0;

    for(int i=0; i<6*num_star; i++)
    {
        if(conditions[i % 6] == 1)
            result_vector[tmp++] = three_stars[i];
    }
    result_vector[size*num_star] = 8e36;

    gauss_newton(result_vector, num_star, conditions);
}
