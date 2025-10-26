#include <cmath>

#include "diff.hpp"
#include "constans.hpp"
#include "transform.hpp"

#define MAX_ITER_NEWTON 100

double solve_kepler_eq(double e, double M) {
    // Solve Kepler�s Equation for the eccentric anomaly
    return newtons_method(M, e, M);
}

void derivative_by_m(double dt, kepler_orbit_denorm* denorm, double R_0, double m, double* ra, double* dec) {
    // -- CONST --
    double norm_a = (denorm->a * M_PI * (R_0 / 648000.0)) * LIGHT_YEAR;
    double norm_e = denorm->e;
    double norm_w = denorm->w * M_PI / 180.0;
    double norm_omega = denorm->omega * M_PI / 180.0;
    double norm_i = denorm->i * M_PI / 180.0;
    double delta_t = 365.25 * 86400 * dt;

    // Матричные коэффициенты
    double A0 = cos(norm_w) * cos(norm_omega) - sin(norm_w) * cos(norm_i) * sin(norm_omega);
    double B0 = sin(norm_w) * cos(norm_omega) + cos(norm_w) * cos(norm_i) * sin(norm_omega);
    double A1 = cos(norm_w) * sin(norm_omega) + sin(norm_w) * cos(norm_i) * cos(norm_omega);
    double B1 = cos(norm_w) * cos(norm_i) * cos(norm_omega) - sin(norm_w) * sin(norm_omega);

    // Производная среднего движения по массе
    double u_by_m = 0.5 * sqrt(G / (m * pow(norm_a, 3)));

    double M_by_m = delta_t * u_by_m;

    double M = fmod(sqrt( (m/pow(norm_a, 3)) * G) *delta_t, 2*M_PI);
    double E = solve_kepler_eq(norm_e, M);

    double E_by_m = M_by_m / (1 - norm_e * E);
    double D_by_m = norm_a * norm_e * sin(E) * E_by_m;
    double v_by_m = E_by_m * sqrt(1 - pow(norm_e, 2))/(1 - norm_e * cos(E));

    double v = 2*atan2(
                sqrt(1+norm_e) * sin(E/2),
                sqrt(1-norm_e) * cos(E/2)
            );
    double D = norm_a * (1 - norm_e * cos(E));

    double x_by_m = D_by_m * cos(v) - D * sin(v) * v_by_m;
    double y_by_m = D_by_m * sin(v) + D * cos(v) * v_by_m;

    double pos1_by_m = x_by_m * A0 - y_by_m * B0;
    double pos0_by_m = x_by_m * A1 - y_by_m * B1;

    double pos2 = (double) R_BH_LY * (double) LIGHT_YEAR;

    *ra = (180.0 * 3600.0) / (pos2 * M_PI) * pos1_by_m;
    *dec = (180.0 * 3600.0) / (pos2 * M_PI) * pos0_by_m;
}
