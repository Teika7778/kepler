#include <iostream>
#include "chol.hpp"

// Функция для вычисления индекса в одномерном массиве
int sym_matrix_index(int i, int j, int size) {
    // Только верхний треугольник
    if (i > j) {
        int t = i;
        i = j;
        j = t;
    }
    int sum = 0, k=0;
    while(k != i)
    {
        sum += (size-k);
        k += 1;
    }
    return sum + (j-i);
}

// Получить значение элемента матрицы
double sym_matrix_at(double* matrix, int i, int j, int size) {
    return matrix[sym_matrix_index(i, j, size)];
}

// Изменить значение элемента матрицы
void sym_matrix_change(double* matrix, int i, int j, int size, double value) {
    matrix[sym_matrix_index(i, j, size)] = value;
}


void find_inverse(double** A, double** res, int size)
{
    for (int i=0; i<size; i++)
    {
        double answ[size];
        double b[size];
        for(int j=0; j<size; j++)
        {
            if (j == i) b[j] = 1;
            else b[j] =0;
        }
        solve_eq(A, b, size, answ);

        for (int j=0; j<size; j++)
        {
            res[j][i] = answ[j];
        }
    }
}