#pragma once

#include "struct.hpp"

void derivative_by_m(double dt, kepler_orbit_denorm* denorm, double R_0, double m, double* ra, double* dec);
