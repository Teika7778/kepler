#include <iostream>

// Функция для вычисления индекса в одномерном массиве
int sym_matrix_index(int i, int j, int size) {
    // Только верхний треугольник
    if (i > j) {
        int t = i;
        i = j;
        j = t;
    }
    // Формула для индекса в нижнем треугольнике (включая диагональ)
    return i * (i + 1) / 2 + j;
}

// Получить значение элемента матрицы
double sym_matrix_at(double* matrix, int i, int j, int size) {
    return matrix[sym_matrix_index(i, j, size)];
}

// Изменить значение элемента матрицы
void sym_matrix_change(double* matrix, int i, int j, int size, double value) {
    matrix[sym_matrix_index(i, j, size)] = value;
}


int main()
{
    double matrix[9];

    for (int i=0; i<9; i++) matrix[i] = i;

    std::cout << matrix[2*3 +1] << std::endl;
}