#include <cmath>
#include "chol.hpp"

void decompose(double** A, double** L, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            double sum = 0.0;
            for (int k = 0; k < j; k++) {
                sum += L[i][k] * L[j][k];
            }
            if (i == j) {
                L[i][j] = sqrt(A[i][i] - sum);
            } else {
                L[i][j] = (A[i][j] - sum) / L[j][j];
            }
        }
    }
}

void solve_eq(double** A, double* b, int n, double* x) {
    double** L = new double*[n];
    for (int i = 0; i < n; i++) {
        L[i] = new double[n]();
    }
    decompose(A, (double**) L, n);
    // ����� ��������� L LT x = b
    // ���������� ��� ��� L y = b, ��� y = LT x
    double* y = new double[n];
    for (int i = 0; i < n; i++) {
        double s = 0;
        for (int j = 0; j < i; j ++) {
            s += L[i][j] * y[j];
        }
        y[i] = (b[i] - s) / L[i][i];
    }
    // ������ ����� ��������� LT x = y
    for (int i = n - 1; i >= 0; i--) {
        double s = 0.0;
        for (int j = i + 1; j < n; j++) {
            s += L[j][i] * x[j]; // L^T[i][j] = L[j][i]
        }
        x[i] = (y[i] - s) / L[i][i];
    }
    for (int i = 0; i < n; i++) {
        delete[] L[i];
    }
    delete[] L;
    delete[] y;
}
