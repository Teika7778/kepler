#include <stdlib.h>
#include <cmath>

#include "integration.hpp"
#include "transform.hpp"
#include "constans.hpp"
#include "normalize.hpp"

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
    double d = (double) R_BH_LY * (double) LIGHT_YEAR;
    double c = 180 / M_PI * 3600;
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
            x[i * STATE_SIZE + j] = pos[j];
            x[i * STATE_SIZE + j + 3] = velo[j];
        }
    }
    x[0] *= c/d;
    x[1] *= c/d;
    x[6] *= c/d;
    x[7] *= c/d;
    x[12] *= c/d;
    x[13] *= c/d;


}

void dxdt(double t, double* x, double* xdot, void* data)
{
    double d = (double) R_BH_LY * (double) LIGHT_YEAR;
    double c = 180 / M_PI * 3600;
    struct simulation_data* sim_data = (struct simulation_data*)(data);

    for (int i=0; i<sim_data->NBODIES; i++)
    {
        xdot[i*STATE_SIZE] = c/d*x[i*STATE_SIZE+3];
        xdot[i*STATE_SIZE+1] = c/d*x[i*STATE_SIZE+4];
        xdot[i*STATE_SIZE+2] = x[i*STATE_SIZE+5];


        double rx = x[i*STATE_SIZE] * d / c  ;  // относительно центра (черной дыры)
        double ry = x[i*STATE_SIZE+1] * d / c ;
        double rz = x[i*STATE_SIZE+2];
        //double r3 = pow(d, 3);
        double r3 = pow(sqrt(rx*rx+ry*ry+rz*rz), 3);
        xdot[i*STATE_SIZE+3] = (-sim_data->Grav*sim_data->M_bh * rx) / r3;
        xdot[i*STATE_SIZE+4] = (-sim_data->Grav*sim_data->M_bh * ry) / r3;
        xdot[i*STATE_SIZE+5] = (-sim_data->Grav*sim_data->M_bh * rz) / r3;
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
