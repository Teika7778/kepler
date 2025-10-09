#pragma once

#include "struct.hpp"

// funciton takes kepler denoemalized orbit and transforms data to normalized form
// Takes denorm with data, norm without data. R_0 is expected to be in Light Years, M_0 in kg
void normalize(kepler_orbit_denorm* denorm, kepler_orbit* norm, double R_0, double M_0);

