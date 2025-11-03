#include <cmath>
#include <iostream>

#include "constans.hpp"
#include "struct.hpp"
#include "normalize.hpp"
#include "transform.hpp"
#include "diff.hpp"
#include "gauss_newton.hpp"
#include "integration.hpp"


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

    // init states

    FILE* file_s2 = fopen("integration_s2.txt", "w");
    FILE* file_s38 = fopen("integration_s38.txt", "w");
    FILE* file_s55 = fopen("integration_s55.txt", "w");

    init_states(x);
    init_deriv(dxdm);

    double dt = 300000;   // Шаг симуляции

    struct simulation_data_star data = {G, M_BH, NBODIES};  // Дополнительные данные для ode

    double* x_for_deriv = (double*)malloc(sizeof(double)*NBODIES);
    double* y_for_deriv = (double*)malloc(sizeof(double)*NBODIES);
    double* z_for_deriv = (double*)malloc(sizeof(double)*NBODIES);

    struct simulation_data_deriv data_deriv = {G, M_BH, NBODIES, x_for_deriv, y_for_deriv, z_for_deriv};

    while (simulationTime < 1e9) // изменить
    {
        //std::cout << x[0] << " " << x[1] << " " << x[2] <<std::endl;
        fprintf(file_s2, "%lf %lf\n", x[0], x[1]);
        fprintf(file_s38, "%lf %lf\n", x[6], x[7]);
        fprintf(file_s55, "%lf %lf\n", x[12], x[13]);
        ode(&rk_4, x, NBODIES*STATE_SIZE_STAR, simulationTime, simulationTime+dt, dxdt, &data);
        for (int i=0; i < NBODIES; i++)
        {
            data_deriv.x[i] = x[i*STATE_SIZE_STAR];
            data_deriv.y[i] = x[i*STATE_SIZE_STAR+1];
            data_deriv.z[i] = x[i*STATE_SIZE_STAR+2];
        }
        ode(&rk_4, dxdm, NBODIES*STATE_SIZE_DERIV, simulationTime, simulationTime+dt, dxdmdt, &data_deriv);
        printf("%.10e %.10e %.10e\n", dxdm[0], dxdm[1], dxdm[2]);
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
}
