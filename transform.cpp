#include "transform.h"
#include <iostream>
#include <cmath>
#define MAX_ITER_NEWTON 100

// Newton
double f(double E, double e, double M)
{
    // Решаемое уравнение
    return E - e*sin(E) - M;
}

double f_der(double E, double e)
{
    // Производная по E
    return 1 - e*cos(E);
}
double newtons_method(double E_0, double e, double M)
{
    double En = E_0;
    int i = 0;

    double arr[MAX_ITER_NEWTON];

    arr[0] = En;

    while (++i != MAX_ITER_NEWTON)
    {
        arr[i] = arr[i-1] - f(arr[i-1], e, M) / f_der(arr[i-1], e);

        if (arr[i] - arr[i-1] == 0)
            return arr[i];
        if (f(arr[i], e, M) == 0)
            return arr[i];


        if (std::isnan(arr[i]) || std::isinf(arr[i])) return arr[i];

        for (int j=0; j<i; j++)
        {
            if (arr[j] == arr[i])
            {
                double best_root = arr[j];
                for (int m = j; m<i; m++)
                {
                    if (fabs(f(arr[m], e, M)) < fabs(f(best_root, e, M)))
                        best_root =  arr[m];
                }

                return best_root;
            }
        }
    }

    return arr[MAX_ITER_NEWTON-1];
}

// Newton END

double calc_M(kepler_orbit* orbit, double grav_param, double t) {
    // Determine the time difference dt in seconds with
    double delta_t = 86400 * (t - orbit->t0);
    // Calculate mean anomaly M(t) from
    double M = orbit->M0 + delta_t * sqrt(grav_param / pow(orbit->a, 3));
    // Normalize M(t) to be in [0; 2pi)
    return fmod(M, 2*M_PI);
}

double solve_kepler_eq(kepler_orbit* orbit, double M) {
    // Solve Kepler’s Equation for the eccentric anomaly
    return newtons_method(M, orbit->e, M);
}

double true_anomaly(kepler_orbit* orbit, double E) {
    double e = orbit->e;
    return 2*atan2(
                sqrt(1+e) * sin(E/2),
                sqrt(1-e) * cos(E/2)
            );
}

double calc_dist(kepler_orbit* orbit, double E) {
    return orbit->a * (1 - orbit->e * cos(E));
}

void calc_pos(double dist, double true_anomaly, double* ret) {
    ret[0] = dist * cos(true_anomaly);
    ret[1] = dist * sin(true_anomaly);
    ret[2] = 0;
}

void calc_velo(kepler_orbit* orbit, double dist, double E, double grav_param, double* ret) {
    ret[0] = (sqrt(grav_param * orbit->a) / dist) * (-1 * sin(E));
    ret[1] = (sqrt(grav_param * orbit->a) / dist) * (sqrt(1-pow(orbit->e,2)) * cos(E));
    ret[2] = 0;
}

void transform_cords(double* s, kepler_orbit* orbit, double* ret) {
    double w = orbit->w;
    double i = orbit->i;
    double O = orbit->omega;
    double x = s[0]; double y = s[1];
    ret[0] =  x*(cos(w)*cos(O)-sin(w)*cos(i)*sin(O))
            - y*(sin(w)*cos(O)+cos(w)*cos(i)*sin(O));

    ret[1] =  x*(cos(w)*sin(O)+sin(w)*cos(i)*cos(O))
            + y*(cos(w)*cos(i)*cos(O)-sin(w)*sin(O));

    ret[2] = x*sin(w)*sin(i) + y*cos(w)*sin(i);
}

void kepler_to_cart(kepler_orbit* orbit, double t, double grav_param,
                    double* pos, double* velo) {

    double M = calc_M(orbit, grav_param, t);
    double E = solve_kepler_eq(orbit, M);
    double v = true_anomaly(orbit, E);
    double d = calc_dist(orbit, E);
    double r[3];
    double u[3];
    calc_pos(d, v, r);
    calc_velo(orbit, d, E, grav_param, u);
    transform_cords(r, orbit, pos);
    transform_cords(u, orbit, velo);

}
