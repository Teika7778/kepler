#pragma once

int sym_matrix_index(int i, int j, int size);

double sym_matrix_at(double* matrix, int i, int j, int size);

void sym_matrix_change(double* matrix, int i, int j, int size, double value);