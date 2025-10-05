#include "transform.h"
#include "normalize.h"
#define G 6.67428e-11
#include <cmath>

int main() {
    kepler_orbit_denorm denorm;

    denorm.a = 0.121;
    denorm.e = 0.872;
    denorm.w = 68.9;
    denorm.omega = 231.9;
    denorm.i = 138.1;
    denorm.T0 = 2002.27;
    denorm.t0 = 2002.9;

    double R_0 = 9460730472580800.0 * 26996; // distance in meters
    double M_0 = 8.54e36; // mass in kg

    kepler_orbit n;

    int t0 = 5;
    int step = 5000;
    int steps = 6;
    int m_val = t0 + step*steps;

    double velo[3];
    double pos[3];
    double ra, dec;

    for (int t = t0; t < m_val; t += step) {
        denorm.t0 = 2002.9 + (double) t;
        normalize(&denorm, &n, R_0, M_0);
        double grav = sqrt(M_0*G/pow(n.a, 3));
        kepler_to_cart(&n, t, grav, pos, velo);


        ra = atan2(pos[2], sqrt(pow(pos[0],2) + pow(pos[1],2)));
        dec = atan2(pos[1], pos[0]);
        std::cout << ra << " " << dec << std::endl;
    }
}
