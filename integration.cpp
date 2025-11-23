#include <stdlib.h>
#include <cmath>
#include <iostream>

#include "integration.hpp"
#include "transform.hpp"
#include "constans.hpp"
#include "normalize.hpp"



void dxdt(double t, double* x, double* xdot, void* data)
{
    struct simulation_data_star* sim_data = (struct simulation_data_star*)(data);

    for (int i=0; i<sim_data->NBODIES; i++)
    {
        xdot[i*STATE_SIZE_STAR]   = x[i*STATE_SIZE_STAR+3];
        xdot[i*STATE_SIZE_STAR+1] = x[i*STATE_SIZE_STAR+4];
        xdot[i*STATE_SIZE_STAR+2] = x[i*STATE_SIZE_STAR+5];


        double rx = x[i*STATE_SIZE_STAR];  // относительно центра (черной дыры)
        double ry = x[i*STATE_SIZE_STAR+1];
        double rz = x[i*STATE_SIZE_STAR+2];

        double r3 = pow(sqrt(rx*rx + ry*ry + rz*rz), 3);

        xdot[i*STATE_SIZE_STAR+3] = (-sim_data->Grav*sim_data->M_bh * rx) / r3;
        xdot[i*STATE_SIZE_STAR+4] = (-sim_data->Grav*sim_data->M_bh * ry) / r3;
        xdot[i*STATE_SIZE_STAR+5] = (-sim_data->Grav*sim_data->M_bh * rz) / r3;
    }
}


void ode(rk4* self, double* x, int n, double t0, double t1, void (*f)(double, double*, double*, void*), void* data)
{

    if (self->k1 == NULL)
        self->k1 = (double*)malloc(sizeof(double)*n);

    if (self->k2 == NULL)
        self->k2 = (double*)malloc(sizeof(double)*n);

    if (self->k3 == NULL)
        self->k3 = (double*)malloc(sizeof(double)*n);

    if (self->k4 == NULL)
        self->k4 = (double*)malloc(sizeof(double)*n);

    if (self->tmp == NULL)
        self->tmp = (double*)malloc(sizeof(double)*n);

    double h = t1-t0;

    f(t0, x, self->k1, data); // now k1 has the right part of f(t, x)

    for (int i=0; i<n; i++)
        self->tmp[i] = x[i] + h*0.5*self->k1[i]; // now tmp is x_n + h/2*k1

    f(t0+h*0.5, self->tmp, self->k2, data);  // now k2 has the right part f(t/2, x_n+h/2*k1)

    for (int i=0; i<n; i++)
        self->tmp[i] = x[i] + h*0.5*self->k2[i];  // now tmp is x_n + h/2*k2

    f(t0+h*0.5, self->tmp, self->k3, data);  // now k3 has the right part f(t/2, x_n+h/2*k2)

    for (int i=0; i<n; i++)
        self->tmp[i] = x[i] + h*self->k3[i]; // now tmp is x_n + h*k3

    f(t0, self->tmp, self->k4, data);  // now fx has the right part f(k3)

    for (int i=0; i<n; i++)
        x[i] +=  1./6.*h*(self->k1[i]+2*self->k2[i]+2*self->k3[i]+self->k4[i]);

}

void rk4Free(rk4* rk)
{
    free(rk->k1);
    free(rk->k2);
    free(rk->k3);
    free(rk->k4);
    free(rk->tmp);
}


void wrap_integration(double* x, double t, double M_bh, rk4 rk_4)
{

    double dt = 86400;   // Шаг - день

    if (t<0){
        dt *= -1;
    }

    struct simulation_data_star data = {G, M_bh, 1};  // Дополнительные данные для ode

    double local_time = 0;

    while (std::abs(local_time) <= std::abs(t))
    {
        ode(&rk_4, x, 1*STATE_SIZE_STAR, local_time, local_time+dt, dxdt, &data);
        local_time += dt;  // Увеличиваем время симуляции
    }

}
