#pragma once

double gauss_newton(kepler_orbit_denorm* stars, double M_bh, int steps=20, bool numerical=true);