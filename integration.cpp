#include <stdlib.h>
#include <cmath>
#include <iostream>

#include "integration.hpp"
#include "transform.hpp"
#include "constans.hpp"
#include "normalize.hpp"

void init_star_state(double* x, kepler_orbit_denorm orbit_denorm, double M_bh)
{
    double pos[3];
    double velo[3];
    kepler_orbit orbit;
    orbit_denorm.t0 = orbit_denorm.T0 ;
    normalize(&orbit_denorm, &orbit, R_BH_LY, M_bh);
    double grav = M_bh * G;
    kepler_to_cart(&orbit, grav, pos, velo);
    for (int j = 0; j < 3; j++)
    {
        x[j] = pos[j];
        x[j + 3] = velo[j];
    }

    for (int j=6; j<12; j++)
        x[j] = 0; // init deriv
}


void dxdt(double* x, double* xdot, void* data)
{
    simulation_data_star* sim_data = (simulation_data_star*)data;

    // Звезда
    xdot[0] = x[3];
    xdot[1] = x[4];
    xdot[2] = x[5];

    double rx = x[0];  // относительно центра (черной дыры)
    double ry = x[1];
    double rz = x[2];
    double r3 = pow(sqrt(rx*rx + ry*ry + rz*rz), 3);
    
    xdot[3] = (-sim_data->Grav*sim_data->M_bh * rx) / r3;
    xdot[4] = (-sim_data->Grav*sim_data->M_bh * ry) / r3;
    xdot[5] = (-sim_data->Grav*sim_data->M_bh * rz) / r3;

    //Производная
    xdot[6] = x[9];
    xdot[7] = x[10];
    xdot[8] = x[11];

    double x_sim = x[0];
    double y_sim = x[1];
    double z_sim = x[2];

    double dxdm = x[6];
    double dydm = x[7];
    double dzdm = x[8];

    double r = sqrt(pow(x_sim, 2) + pow(y_sim, 2) + pow(z_sim, 2));

    r3 = pow(r, 3);
    double r5 = pow(r, 5);

    // Коэффициенты матрицы Якоби
    double Axx = -sim_data->Grav * sim_data->M_bh * (1.0/r3 - 3.0*x_sim*x_sim/r5);
    double Axy = 3.0 * sim_data->Grav * sim_data->M_bh * x_sim * y_sim / r5;
    double Axz = 3.0 * sim_data->Grav * sim_data->M_bh * x_sim * z_sim / r5;

    double Ayx = Axy;  // Симметрия
    double Ayy = -sim_data->Grav * sim_data->M_bh * (1.0/r3 - 3.0*y_sim*y_sim/r5);
    double Ayz = 3.0 * sim_data->Grav * sim_data->M_bh * y_sim * z_sim / r5;

    double Azx = Axz;  // Симметрия
    double Azy = Ayz;  // Симметрия
    double Azz = -sim_data->Grav * sim_data->M_bh * (1.0/r3 - 3.0*z_sim*z_sim/r5);

    // Прямые производные силы по массе
    double dFxdM = -sim_data->Grav * x_sim / r3;
    double dFydM = -sim_data->Grav * y_sim / r3;
    double dFzdM = -sim_data->Grav * z_sim / r3;
    // Уравнения для производных скоростей по массе
    // d(dvx/dM)/dt = Axx*dxdm + Axy*dydm + Axz*dzdm + dFxdM
    xdot[9] =   Axx * dxdm +  Axy * dydm + Axz * dzdm + dFxdM;
    // d(dvy/dM)/dt = Ayx*dxdm + Ayy*dydm + Ayz*dzdm + dFydM
    xdot[10] =  Ayx * dxdm + Ayy * dydm + Ayz * dzdm + dFydM;
    // d(dvz/dM)/dt = Azx*dxdm + Azy*dydm + Azz*dzdm + dFzdM
    xdot[11] =  Azx * dxdm + Azy * dydm + Azz * dzdm + dFzdM;
    
}


void ode(rk4* self, double* x, int n, double t0, double t1, void (*f)(double*, double*, void*), void* data)
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

    f(x, self->k1, data); // now k1 has the right part of f(t, x)

    for (int i=0; i<n; i++)
        self->tmp[i] = x[i] + h*0.5*self->k1[i]; // now tmp is x_n + h/2*k1

    f(self->tmp, self->k2, data);  // now k2 has the right part f(t/2, x_n+h/2*k1)

    for (int i=0; i<n; i++)
        self->tmp[i] = x[i] + h*0.5*self->k2[i];  // now tmp is x_n + h/2*k2

    f(self->tmp, self->k3, data);  // now k3 has the right part f(t/2, x_n+h/2*k2)

    for (int i=0; i<n; i++)
        self->tmp[i] = x[i] + h*self->k3[i]; // now tmp is x_n + h*k3

    f(self->tmp, self->k4, data);  // now fx has the right part f(k3)

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

    struct simulation_data_star data = {G, M_bh};  // Дополнительные данные для ode

    double local_time = 0;

    while (std::abs(local_time) <= std::abs(t))
    {
        ode(&rk_4, x, STATE_SIZE_STAR, local_time, local_time+dt, dxdt, &data);
        local_time += dt;  // Увеличиваем время симуляции
    }

}
