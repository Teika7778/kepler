#include <cmath>
#include <iostream>

#include "constans.hpp"
#include "struct.hpp"
#include "normalize.hpp"
#include "transform.hpp"

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


    int fractioning = 3000;

    double delta_t = 30. / fractioning;

    kepler_orbit orbit;

    double velo[3];
    double pos[3];
    double cur_ra, cur_dec;

    double d = R_BH_LY;

    double BH_x = (d * cos(DEC_BH) * cos(RA_BH) ) * LIGHT_YEAR;
    double BH_y = (d * cos(DEC_BH) * sin(RA_BH) ) * LIGHT_YEAR;
    double BH_z = (d * sin(DEC_BH) )* LIGHT_YEAR ;

    FILE* files[3];
    files[0] = fopen("s2_angles.txt", "w");
    files[1] = fopen("s38_angles.txt", "w");
    files[2] = fopen("s55_angles.txt", "w");

    if (!files[0] || !files[1] || !files[2]) return 0;

    for (size_t star=0; star<1; star++)
    {

        for (size_t i = 0; i< fractioning; i++)
        {

            stars_denorm[star].t0 = stars_denorm[star].T0 + i*delta_t;
            normalize(&stars_denorm[star], &orbit, R_BH_LY, M_BH);

            // double grav = sqrt( (M_BH/pow(orbit.a, 3)) * G );
            double grav = M_BH * G;

            kepler_to_cart(&orbit, grav, pos, velo);

            pos[0] += 0;
            pos[1] += 0;
            pos[2] += d * LIGHT_YEAR; // Можно опустить

            // cur_dec = atan2(pos[2], sqrt(pow(pos[0],2) + pow(pos[1],2))) * 180.0 / M_PI;
            // cur_ra = atan2(pos[1], pos[0]) * 180.0 / M_PI;
            cur_dec = pos[0]/pos[2] * 180.0 / M_PI;
            cur_ra = pos[1]/pos[2] * 180.0 / M_PI;

            // fprintf(files[star], "%.15f %.15f\n", (cur_ra + 360 - RA_BH* 180.0 / M_PI) * 3600, (cur_dec - DEC_BH * 180.0 / M_PI) * 3600);
            fprintf(files[star], "%.15f %.15f\n", cur_ra*3600, cur_dec*3600);
            printf("%.15f %.15f %.15f\n",stars_denorm[star].T0 + i*delta_t ,cur_ra*3600, cur_dec*3600);

        }

    }

    fclose(files[0]);
    fclose(files[1]);
    fclose(files[2]);

    return 0;
}
