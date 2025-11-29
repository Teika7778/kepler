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
    double parameters[6];

    double EPS = 5;

    // Хард код, параметры s2 на первом наблюдении
    parameters[0] = -1.34542432318288867188e+13 - -1.34542432318288867188e+13/ EPS;
    parameters[1] = 2.74369376973100244141e+12 + 2.74369376973100244141e+12 / EPS;
    parameters[2] = 1.17902665512916289062e+13 - 1.17902665512916289062e+13 / EPS;
    parameters[3] = 9.63030905550202805898e+03 + 9.63030905550202805898e+03 / EPS;
    parameters[4] = 2.38743901750857585284e+04 - 2.38743901750857585284e+04 / EPS;
    parameters[5] = 1e32;


    kepler_orbit_denorm denorm_orbit_s2 =
    {
        0.126,  // a
        0.884,  // e
        71.36,   // w
        234.50,  // omega
        136.78,  // i
        2002.32, // T0
        2002.578  //t0
    };
    kepler_orbit_denorm denorm_orbit_s38 =
    {
        0.140,  // a
        0.818,  // e
        18.4,   // w
        101.8,  // omega
        166.22,  // i
        2003.30, // T0
        2003.30  //t0
    };
    kepler_orbit_denorm denorm_orbit_s55 =
    {
        0.109,  // a
        0.74,  // e
        133.5,   // w
        129.9,  // omega
        141.7,  // i
        2009.31, // T0
        2009.31  //t0
    };
    kepler_orbit_denorm stars_denorm[3] =
    {
        denorm_orbit_s2,
        denorm_orbit_s38,
        denorm_orbit_s55
    };

    kepler_orbit orbit;

    normalize(&denorm_orbit_s2, &orbit, (double)R_BH_LY*LIGHT_YEAR, M_BH);

    double pos[3];
    double velo[3];

    kepler_to_cart(&orbit, G, pos, velo);

    std::cout << pos[0] <<  " " << pos[1] <<  " " << pos[2] << std::endl;
    std::cout << velo[0] <<  " " << velo[1] <<  " " << velo[2];

    //gauss_newton(parameters);
}