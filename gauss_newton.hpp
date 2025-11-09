#pragma once

#define MAX_ITER_GAUSS_NEWTON 11

double gauss_newton(kepler_orbit_denorm* stars, double M_bh);

double gauss_newton_2(kepler_orbit_denorm* stars, double M_bh);