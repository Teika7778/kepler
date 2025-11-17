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

    double d = (double) R_BH_LY * (double) LIGHT_YEAR;
    double c = 180. / M_PI * 3600.;

    x[0] *= c/d;
    x[1] *= c/d;
}


void init_states(double* x)
{

        kepler_orbit_denorm denorm_orbit_s2 =
    {
        0.126,  // a
        0.884,  // e
        71.36,   // w
        234.50,  // omega
        136.78,  // i
        2002.32, // T0
        2002.32  //t0
    };
    kepler_orbit_denorm denorm_orbit_s38 =
    {
        0.140,  // a
        0.818,  // e
        18.4,   // w
        101.8,  // omega
        166.22,  // i
        2003.30, // T0
        2003.30  //t0
    };
    kepler_orbit_denorm denorm_orbit_s55 =
    {
        0.109,  // a
        0.74,  // e
        133.5,   // w
        129.9,  // omega
        141.7,  // i
        2009.31, // T0
        2009.31  //t0
    };
    kepler_orbit_denorm stars_denorm[3] =
    {
        denorm_orbit_s2,
        denorm_orbit_s38,
        denorm_orbit_s55
    };

    for(int i=0; i<3; i++)
    {

        double pos[3];
        double velo[3];
        kepler_orbit orbit;

        stars_denorm[i].t0 = stars_denorm[i].T0 ;
        normalize(&stars_denorm[i], &orbit, R_BH_LY, M_BH);

        double grav = M_BH * G;

        kepler_to_cart(&orbit, grav, pos, velo);

        for (int j = 0; j < 3; j++)
        {
            x[i * STATE_SIZE_STAR + j] = pos[j];
            x[i * STATE_SIZE_STAR + j + 3] = velo[j];
        }
    }

    double d = (double) R_BH_LY * (double) LIGHT_YEAR;
    double c = 180. / M_PI * 3600.;

    x[0] *= c/d;
    x[1] *= c/d;
    x[6] *= c/d;
    x[7] *= c/d;
    x[12] *= c/d;
    x[13] *= c/d;
}

void dxdmdt(double t, double* x, double* xdot, void* data)
{
    struct simulation_data_deriv* sim_data = (struct simulation_data_deriv*)(data);

    double d = (double) R_BH_LY * (double) LIGHT_YEAR;
    double c = 180. / M_PI * 3600.;

    for (int i=0; i<sim_data->NBODIES; i++)
    {
        xdot[i*STATE_SIZE_DERIV]   = x[i*STATE_SIZE_DERIV+3];
        xdot[i*STATE_SIZE_DERIV+1] = x[i*STATE_SIZE_DERIV+4];
        xdot[i*STATE_SIZE_DERIV+2] = x[i*STATE_SIZE_DERIV+5];

        double x_sim = d/c * sim_data->x[i];
        double y_sim = d/c * sim_data->y[i];
        double z_sim = sim_data->z[i];

        double dxdm = x[i*STATE_SIZE_DERIV];
        double dydm = x[i*STATE_SIZE_DERIV+1];
        double dzdm = x[i*STATE_SIZE_DERIV+2];

        double r = sqrt(pow(x_sim, 2) + pow(y_sim, 2) + pow(z_sim, 2));

        double r3 = pow(r, 3);
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
        xdot[i*STATE_SIZE_DERIV+3] =   Axx * dxdm +  Axy * dydm + Axz * dzdm + dFxdM;

        // d(dvy/dM)/dt = Ayx*dxdm + Ayy*dydm + Ayz*dzdm + dFydM
        xdot[i*STATE_SIZE_DERIV+4] =  Ayx * dxdm + Ayy * dydm + Ayz * dzdm + dFydM;

        // d(dvz/dM)/dt = Azx*dxdm + Azy*dydm + Azz*dzdm + dFzdM
        xdot[i*STATE_SIZE_DERIV+5] =  Azx * dxdm + Azy * dydm + Azz * dzdm + dFzdM;
    }
}


void dxdt(double t, double* x, double* xdot, void* data)
{
    struct simulation_data_star* sim_data = (struct simulation_data_star*)(data);

    double d = (double) R_BH_LY * (double) LIGHT_YEAR;
    double c = 180 / M_PI * 3600;

    for (int i=0; i<sim_data->NBODIES; i++)
    {
        xdot[i*STATE_SIZE_STAR] = c/d * x[i*STATE_SIZE_STAR+3];
        xdot[i*STATE_SIZE_STAR+1] = c/d * x[i*STATE_SIZE_STAR+4];
        xdot[i*STATE_SIZE_STAR+2] = x[i*STATE_SIZE_STAR+5];


        double rx = x[i*STATE_SIZE_STAR] * d /c;  // относительно центра (черной дыры)
        double ry = x[i*STATE_SIZE_STAR+1] * d/ c;
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


void wrap_integration(double* x, double* deriv, double t, double M_bh, rk4 rk_4, double** arrays_for_deriv)
{

    double dt = 86400;   // Шаг - день

    if (t<0){
        dt *= -1;
    }

    struct simulation_data_star data = {G, M_bh, 1};  // Дополнительные данные для ode

    double* x_for_deriv = arrays_for_deriv[0];
    double* y_for_deriv = arrays_for_deriv[1];
    double* z_for_deriv = arrays_for_deriv[2];

    struct simulation_data_deriv data_deriv = {G, M_bh, 1, x_for_deriv, y_for_deriv, z_for_deriv};

    double local_time = 0;

    while (std::abs(local_time) <= std::abs(t))
    {
        ode(&rk_4, x, 1*STATE_SIZE_STAR, local_time, local_time+dt, dxdt, &data);
        data_deriv.x[0] = x[0];
        data_deriv.y[0] = x[1];
        data_deriv.z[0] = x[2];
        ode(&rk_4, deriv, 1*STATE_SIZE_DERIV, local_time, local_time+dt, dxdmdt, &data_deriv);
        local_time += dt;  // Увеличиваем время симуляции
    }

}
