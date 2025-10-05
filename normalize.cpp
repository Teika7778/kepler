#include "normalize.h"
#include "transform.h"
#define G 6.67428e-11


void normalize(kepler_orbit_denorm* denorm, kepler_orbit* norm, double R_0, double M_0)
{

    // Normalize Semi-major axis
    // arcsec --> m

    // Если R_0 не в метрах
    // R_0 = 9460730472580800 * R_0 астр.года

    norm->a = denorm->a * M_PI / 648000.0 * R_0;

    // No need to normalize Eccentricity
    norm->e = denorm->e;

    // Normalize Argument of periapsis. (deg -> rad)
    norm->w = denorm->w * M_PI / 180.0;

    // Normalize Longitude of ascending node  (deg -> rad)
    norm->omega = denorm->omega * M_PI / 180.0;

    // Normalize Inclination  (deg -> rad)
    norm->i = denorm->i * M_PI / 180.0;

    // Calculate mean_movement
    double mean_movement = sqrt(M_0*G/pow(norm->a, 3));

    double delta_t = 8400*(denorm->t0 - denorm->T0);
    double mean_anomaly = mean_movement*delta_t;

    //mean_anomaly = fmod(M_0, 2*M_PI); // Возможное место ошибок

    // Normalize Mean anomaly
    norm->M0 = mean_anomaly;

    // No need to normalize epoch
    norm->t0 = denorm->t0;

}