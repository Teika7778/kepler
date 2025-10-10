#include <cmath>
#include <iostream>

#include "constans.hpp"
#include "struct.hpp"
#include "normalize.hpp"
#include "transform.hpp"

// Сделать на эллипсе 1000 точек так, чтобы они полностью покрывали период
// Перевести координаты в декартовы
// Добавить к координатам декартовы координаты Черной дыры
// Сделать 999 шагов и посчитать изменение ra и dec


int main() {
    kepler_orbit_denorm denorm_orbit =
    {
        0.121,  // a
        0.872,  // e
        68.9,   // w
        231.9,  // omega
        138.1,  // i
        2002.27, // T0
        2002.27  //t0
    };

    int fractioning = 1000;

    double delta_t = 16. / fractioning;

    double date = 2002.27;

    kepler_orbit orbit;

    double velo[3];
    double pos[3];
    double cur_ra, cur_dec, prev_ra, prev_dec;

    double d = R_BH_LY * LIGHT_YEAR;
    double BH_x = d * cos(DEC_BH) * cos(RA_BH);
    double BH_y = d * cos(DEC_BH) * sin(RA_BH);
    double BH_z = d * sin(DEC_BH);

    FILE* file = fopen("angles.txt", "w");

    if (!file) return 0;

    for (size_t i = 0; i< fractioning; i++)
    {

        denorm_orbit.t0 = denorm_orbit.T0 + i*delta_t;
        normalize(&denorm_orbit, &orbit, R_BH_LY, M_BH);

        double grav = sqrt( (M_BH/pow(orbit.a, 3)) * G );
        kepler_to_cart(&orbit, grav, pos, velo);
        pos[0] += BH_x;
        pos[1] += BH_y;
        pos[2] += BH_z;

        cur_ra = atan2(pos[1], pos[0]) * 180.0 / M_PI;
        cur_dec = atan2(pos[2], sqrt(pow(pos[0],2) + pow(pos[1],2))) * 180.0 / M_PI;





        if ( i != 0)
        {
            fprintf(file, "%f %f\n", cur_ra + 180 - RA_BH* 180.0 / M_PI, -1 * cur_dec - DEC_BH * 180.0 / M_PI);
            //fprintf(file, "%f %f\n", cur_ra - RA_BH * 180.0 / M_PI, cur_dec - DEC_BH * 180.0 / M_PI);
        }

        prev_ra = cur_ra;
        prev_dec = cur_dec;

    }

    fclose(file);

    return 0;
}
