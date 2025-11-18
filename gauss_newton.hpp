#pragma once

#define MAX_ITER_GAUSS_NEWTON 6

double gauss_newton(kepler_orbit_denorm* stars, double M_bh);

void gauss_newton_3(double* parameters);