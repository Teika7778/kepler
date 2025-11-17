#include <cmath>
#include <iostream>

#include "constans.hpp"
#include "struct.hpp"
#include "normalize.hpp"
#include "transform.hpp"
#include "diff.hpp"
#include "gauss_newton.hpp"
#include "integration.hpp"

/*
void diff_m(double delta_t, double mass, kepler_orbit_denorm orbit_d, double* dRA, double* ddec, double eps) {
    double RA_l, DEC_l, RA_r, DEC_r;
    warp(delta_t, mass-eps, orbit_d, &RA_l, &DEC_l);
    warp(delta_t, mass+eps, orbit_d, &RA_r, &DEC_r);

    // *dRA = (RA_r - RA_l) / 2*eps;
    // *ddec = (DEC_r - DEC_l) / 2*eps;
    *dRA = (RA_r - RA_l) / 2*eps;
    *ddec = (DEC_r - DEC_l) / 2*eps ;
}

#define NBODIES 3 //Количество тел в симуляции


double x[NBODIES*STATE_SIZE_STAR]; // Глобальный массив тел

double dxdm[NBODIES*STATE_SIZE_DERIV];

double simulationTime = 0; // Глобальная переменная времени

rk4 rk_4 = {NULL, NULL, NULL, NULL, NULL};  // Глобальная переменная структуры РК-4

int main() {

    double d = (double) R_BH_LY * (double) LIGHT_YEAR;
    double c = 180 / M_PI * 3600;
    // init states

    FILE* file_s2 = fopen("integration_s2.txt", "w");
    FILE* file_s38 = fopen("integration_s38.txt", "w");
    FILE* file_s55 = fopen("integration_s55.txt", "w");

    init_states(x);
    init_deriv(dxdm);

    double dt = 86400;   // Шаг - неделя

    struct simulation_data_star data = {G, M_BH, NBODIES};  // Дополнительные данные для ode

    double* x_for_deriv = (double*)malloc(sizeof(double)*NBODIES);
    double* y_for_deriv = (double*)malloc(sizeof(double)*NBODIES);
    double* z_for_deriv = (double*)malloc(sizeof(double)*NBODIES);

    struct simulation_data_deriv data_deriv = {G, M_BH, NBODIES, x_for_deriv, y_for_deriv, z_for_deriv};

    while (simulationTime < 25.261*365*86400) // изменить
    {
        //std::cout << x[0] << " " << x[1] << " " << x[2] <<std::endl;

        fprintf(file_s2, "%lf %lf\n", x[1], x[0]);
        fprintf(file_s38, "%lf %lf\n", x[7], x[6]);
        fprintf(file_s55, "%lf %lf\n", x[13], x[12]);
        ode(&rk_4, x, NBODIES*STATE_SIZE_STAR, simulationTime, simulationTime+dt, dxdt, &data);
        for (int i=0; i < NBODIES; i++)
        {
            data_deriv.x[i] = x[i*STATE_SIZE_STAR];
            data_deriv.y[i] = x[i*STATE_SIZE_STAR+1];
            data_deriv.z[i] = x[i*STATE_SIZE_STAR+2];
        }
        ode(&rk_4, dxdm, NBODIES*STATE_SIZE_DERIV, simulationTime, simulationTime+dt, dxdmdt, &data_deriv);
        simulationTime += dt;  // Увеличиваем время симуляции
    }

    rk4Free(&rk_4);

    free(x_for_deriv);
    free(y_for_deriv);
    free(z_for_deriv);

    fclose(file_s2);
    fclose(file_s38);
    fclose(file_s55);

    return 0;
}*/


int main()
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

    gauss_newton(stars_denorm, 1e32);
}