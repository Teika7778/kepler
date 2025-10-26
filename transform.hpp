#pragma once

#include <iostream>
#include <cmath>

#include "struct.hpp"

// kepler equation
double kepler_equation(double E, double e, double M);

// kepler equation derivative
double kepler_equation_derivative(double E, double e);

// newtowns method to solve equation
double newtons_method(double E_0, double e, double M);

double solve_kepler_eq(kepler_orbit* orbit, double M);

double true_anomaly(kepler_orbit* orbit, double E);

double calc_dist(kepler_orbit* orbit, double E);

void calc_pos(double dist, double true_anomaly, double* ret);

void calc_velo(kepler_orbit* orbit, double dist, double E, double grav_param, double* ret);

void transform_cords(double* s, kepler_orbit* orbit, double* ret);

void kepler_to_cart(kepler_orbit* orbit, double grav_param,
                    double* pos, double* velo);

void warp(double delta_t, double mass, kepler_orbit_denorm orbit_d, double* RA, double* dec);