#ifndef NORMALIZE
#define NORMALIZE
#include "normalize.h"
#include "transform.h"

// Здесь даны данные из статьи

struct kepler_orbit_denorm {
    double a;     // Semi-major axis (arcsec)
    double e;     // Eccentricity
    double w;     // Argument of periapsis (deg)
    double omega; // Longitude of ascending node (deg)
    double i;     // Inclination (deg)
    double T0;    // Time of periapsis t_0
    double t0;
} typedef kepler_orbit_denorm;

void normalize(kepler_orbit_denorm* denorm, kepler_orbit* norm, double R_0, double M_0);

#endif
