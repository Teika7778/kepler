#pragma once

#define MAX_ITER_GAUSS_NEWTON 100

void gauss_newton(double* parameters, int star_number, int* conditions);

void init_star(double* x, int star_number);

