#pragma once

#define MAX_ITER_GAUSS_NEWTON 850

double gauss_newton(kepler_orbit_denorm* stars, double M_bh);

double gauss_newton_3(double* parameters);