
#ifndef TRANSFORM
#define TRANSFORM

#include <iostream>
#include <cmath>
#define MAX_ITER_NEWTON 100


struct kepler_orbit {
    double a;     // Semi-major axis
    double e;     // Eccentricity
    double w;     // Argument of periapsis
    double omega; // Longitude of ascending node
    double i;     // Inclination
    double M0;    // Mean anomaly at t_0
    double t0;
} typedef kepler_orbit;

// Newton

double f(double E, double e, double M);

double f_der(double E, double e);

double newtons_method(double E_0, double e, double M);
// Newton END

double calc_M(kepler_orbit* orbit, double grav_param, double t);

double solve_kepler_eq(kepler_orbit* orbit, double M);

double true_anomaly(kepler_orbit* orbit, double E);

double calc_dist(kepler_orbit* orbit, double E);

void calc_pos(double dist, double true_anomaly, double* ret);

void calc_velo(kepler_orbit* orbit, double dist, double E, double grav_param, double* ret);

void transform_cords(double* s, kepler_orbit* orbit, double* ret);

void kepler_to_cart(kepler_orbit* orbit, double t, double grav_param,
                    double* pos, double* velo);

#endif
