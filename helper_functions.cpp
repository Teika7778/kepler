#include <iostream>

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
