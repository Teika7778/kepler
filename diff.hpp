#pragma once

#include "struct.hpp"

void derivative_by_m(double dt, kepler_orbit_denorm* denorm, double R_0, double m, double* ra, double* dec);

void count_diff(double dt, kepler_orbit_denorm denorm, double m, double* ra, double* dec, double eps);

void kepler_deriv(double** deriv_mat, kepler_orbit_denorm denorm_orbit, double M_bh);

void full_analytic(double* deriv_vec, double* params, double t_0);
