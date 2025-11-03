#pragma once

#define STATE_SIZE_STAR 6  // Размер вектора состояния системы

#define STATE_SIZE_DERIV 6

struct kepler_orbit_denorm {
    double a;     // Semi-major axis (arcsec)
    double e;     // Eccentricity
    double w;     // Argument of periapsis (deg)
    double omega; // Longitude of ascending node (deg)
    double i;     // Inclination (deg)
    double T0;    // Time of periapsis t_0
    double t0;
} typedef kepler_orbit_denorm;


struct kepler_orbit {
    double a;     // Semi-major axis (m)
    double e;     // Eccentricity (-)
    double w;     // Argument of periapsis (rad)
    double omega; // Longitude of ascending node (rad)
    double i;     // Inclination (rad)
    double M0;    // Mean anomaly at t_0 (rad)
    double t0;    // Epoch (JD)
} typedef kepler_orbit;


struct rk4{
    double* k1;
    double* k2;
    double* k3;
    double* k4;
    double* tmp;
} typedef rk4;


struct simulation_data_star {
    double Grav;
    double M_bh;
    int NBODIES;
};

struct simulation_data_deriv{
    double Grav;
    double M_bh;
    int NBODIES;
    double* x;
    double* y;
    double* z;
};
