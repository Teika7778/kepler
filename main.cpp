#include <cmath>
#include <iostream>

#include "constans.hpp"
#include "struct.hpp"
#include "normalize.hpp"
#include "transform.hpp"


int main() {
    kepler_orbit_denorm denorm_orbit =
    {
        0.121,  // a
        0.872,  // e
        68.9,   // w
        231.9,  // omega
        138.1,  // i
        2002.27, // T0
        2002.9  //t0 
    };

    kepler_orbit orbit;

    normalize(&denorm_orbit, &orbit, R_BH_LY, M_BH);

    double velo[3];
    double pos[3];
    double ra, dec;

    double grav = sqrt( (M_BH/pow(orbit.a, 3)) * G );
    kepler_to_cart(&orbit, grav, pos, velo);

    ra = atan2(pos[2], sqrt(pow(pos[0],2) + pow(pos[1],2)));
    dec = atan2(pos[1], pos[0]);
    std::cout << ra << " " << dec << std::endl;

    return 0;
}
