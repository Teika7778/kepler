#include <cmath>

#include "struct.hpp"
#include "constans.hpp"

void normalize(kepler_orbit_denorm* denorm, kepler_orbit* norm, double R_0, double M_0)
{

    // Normalize Semi-major axis
    // arcsec --> meters

    norm->a = (denorm->a * M_PI  * (R_0 / 648000.0) ) * LIGHT_YEAR;

    // No need to normalize Eccentricity
    norm->e = denorm->e;

    // Normalize Argument of periapsis. (deg -> rad)
    norm->w = denorm->w * M_PI / 180.0;

    // Normalize Longitude of ascending node  (deg -> rad)
    norm->omega = denorm->omega * M_PI / 180.0;

    // Normalize Inclination  (deg -> rad)
    norm->i = denorm->i * M_PI / 180.0;

    // Calculate mean_movement
    double mean_movement = sqrt( (M_0/pow(norm->a, 3)) * G);

    double delta_t = 8400*(denorm->t0 - denorm->T0);
    double mean_anomaly = mean_movement*delta_t;

    mean_anomaly = fmod(mean_anomaly, 2*M_PI); // Возможное место ошибок

    // Normalize Mean anomaly
    norm->M0 = mean_anomaly;

    // No need to normalize epoch
    norm->t0 = denorm->t0;

}