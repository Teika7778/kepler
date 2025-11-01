#pragma once

#define MAX_ITER_GAUSS_NEWTON 12

void eval(kepler_orbit_denorm* stars, double M_bh);

double gauss_newton(kepler_orbit_denorm* stars, double M_bh);