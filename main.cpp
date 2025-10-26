#include <cmath>
#include <iostream>

#include "constans.hpp"
#include "struct.hpp"
#include "normalize.hpp"
#include "transform.hpp"
#include "diff.hpp"
#include "gauss_newton.hpp"


void diff_m(double delta_t, double mass, kepler_orbit_denorm orbit_d, double* dRA, double* ddec, double eps) {
    double RA_l, DEC_l, RA_r, DEC_r;
    warp(delta_t, mass-eps, orbit_d, &RA_l, &DEC_l);
    warp(delta_t, mass+eps, orbit_d, &RA_r, &DEC_r);

    // *dRA = (RA_r - RA_l) / 2*eps;
    // *ddec = (DEC_r - DEC_l) / 2*eps;
    *dRA = (RA_r - RA_l) / 2*eps;
    *ddec = (DEC_r - DEC_l) / 2*eps ;
}


int main() {


    kepler_orbit_denorm denorm_orbit_s2 =
    {
        0.126,  // a
        0.884,  // e
        71.36,   // w
        234.50,  // omega
        136.78,  // i
        2002.32, // T0
        2002.32  //t0
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


    /*
    int fractioning = 3000;

    double delta_t = 30. / fractioning;
    double cur_ra, cur_dec;
    double d = R_BH_LY;

    FILE* files[3];
    files[0] = fopen("s2_angles.txt", "w");
    files[1] = fopen("s38_angles.txt", "w");
    files[2] = fopen("s55_angles.txt", "w");

    if (!files[0] || !files[1] || !files[2]) return 0;

    for (size_t star=0; star<3; star++)
    {

        for (size_t i = 0; i< fractioning; i++)
        {

            warp(i*delta_t, M_BH, stars_denorm[star], &cur_ra, &cur_dec);
            fprintf(files[star], "%.15f %.15f\n", cur_ra, cur_dec);
            derivative_by_m(i*delta_t,  &stars_denorm[star], d*LIGHT_YEAR, M_BH, &cur_ra, &cur_dec);
            //diff_m(i*delta_t, M_BH, stars_denorm[star], &cur_ra, &cur_dec, 1e+26);
            printf("%e %e\n", cur_ra, cur_dec);
        }
    }

    fclose(files[0]);
    fclose(files[1]);
    fclose(files[2]);
    */

    printf("%5.e", gauss_newton(stars_denorm, M_BH));


    return 0;
}
