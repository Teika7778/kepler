#include "transform.h"
#include "normalize.h"
#include <cmath>

int main() {
    kepler_orbit_denorm d;
    kepler_orbit n;

    double grav = 0;

    //normalize(&d, &n);
    int t0 = 5;
    int step = 5;
    int steps = 6;
    int m_val = t0 + step*steps;

    double velo[3];
    double pos[3];
    double ra, dec;

    for (int t = t0; t < m_val; t += step) {
        kepler_to_cart(&n, t, grav, pos, velo);
        ra = atan2(pos[2], sqrt(pow(pos[0],2) + pow(pos[1],2)));
        dec = atan2(pos[1], pos[0]);
        std::cout << ra << " " << dec << std::endl;
    }
}
