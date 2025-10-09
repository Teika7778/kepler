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

    FILE* file = fopen("angles.txt", "w");

    if (!file) return 0;

    for (size_t i = 0; i< fractioning; i++)
    {
        
        denorm_orbit.t0 = denorm_orbit.T0 + i*delta_t;
        normalize(&denorm_orbit, &orbit, R_BH_LY, M_BH);

        double grav = sqrt( (M_BH/pow(orbit.a, 3)) * G );
        kepler_to_cart(&orbit, grav, pos, velo);

        cur_ra = atan2(pos[2], sqrt(pow(pos[0],2) + pow(pos[1],2))) * 180.0 / M_PI;
        cur_dec = atan2(pos[1], pos[0]) * 180.0 / M_PI;

        

        while (cur_ra - prev_ra > 180) cur_ra -= 360;
        while (cur_ra - prev_ra < -180) cur_ra += 360;
        
        // Корректировка склонения  
        while (cur_dec - prev_dec > 180) cur_dec -= 360;
        while (cur_dec - prev_dec < -180) cur_dec += 360;

        if ( i != 0)
        {
            fprintf(file, "%f %f\n", cur_ra-prev_ra, cur_dec-prev_dec);
        }

        prev_ra = cur_ra;
        prev_dec = cur_dec;

    }

    fclose(file);

    return 0;
}
