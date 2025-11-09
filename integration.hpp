#pragma once

#include <iostream>
#include "struct.hpp"

// Установление начальных координат
void init_states(double* x);

void init_star_state(double* x, kepler_orbit_denorm orbit_denorm, double M_bh);

void init_deriv(double* x);

// Взятие производной у массива тел
void dxdt(double t, double* x, double* xdot, void* data);

// Вычисление изозронной производной
void dxdmdt(double t, double* x, double* xdot, void* data);

// Численное решение ОДУ (Рунге-Кутта 4)
void ode(rk4* self, double* x, int n, double t0, double t1, 
void (*f)(double, double*, double*, void*), void* data); 
// Функция правой части ( время, начальный вектор состояния, конечный вектор состояния, данные)

// Освобождение памяти из под структуры rk4
void rk4Free(rk4* rk);

// Обертка над численным интегрированием
// x - вектор состояния звезды, deriv - вектор состояния производной
// t - время, на которое требуется проинтегрировать систему
// rk_4, array_for_deriv - технические переменные
void wrap_integration(double* x, double* deriv, double t, double M_bh, rk4 rk_4, double** arrays_for_deriv);