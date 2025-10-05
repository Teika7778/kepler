#include "normalize.cpp"
#include "transform.cpp"
#include "transform.h"

int main() {
    
    kepler_orbit_denorm* denorm;

    denorm->a = 0.121;
    denorm->e = 0.872;
    denorm->w = 68.9;
    denorm->omega = 231.9;
    denorm->i = 138.1;
    denorm->T0 = 2002.27;
    denorm->t0 = 2002.9;

    double R_0 = 9460730472580800.0 * 26996; // distance in meters
    double M_0 = 8.54e36; // mass in kg

    kepler_orbit* norm;

    normalize(denorm, norm, R_0, M_0);

    return 0;
}
