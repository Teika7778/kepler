#include <cmath>
#include "transform.hpp"
#include "constans.hpp"
#include "normalize.hpp"

#define MAX_ITER_NEWTON 100

double kepler_equation(double E, double e, double M)
{
    return E - e*sin(E) - M;
}

double kepler_equation_derivative(double E, double e)
{
    return 1 - e*cos(E);
}
double newtons_method(double E_0, double e, double M)
{
    // first approximation
    double En = E_0;
    int i = 0;

    double arr[MAX_ITER_NEWTON];

    arr[0] = En;

    while (++i != MAX_ITER_NEWTON)
    {
        // itaration step
        arr[i] = arr[i-1] - kepler_equation(arr[i-1], e, M) / kepler_equation_derivative(arr[i-1], e);

        if (arr[i] - arr[i-1] == 0) // if value didnt change - no sense going more
            return arr[i];
        if (kepler_equation(arr[i], e, M) == 0) // root found
            return arr[i];

        if (std::isnan(arr[i]) || std::isinf(arr[i])) return arr[i]; // errors found

        for (int j=0; j<i; j++)
        {
            if (arr[j] == arr[i]) // loop found
            {
                double best_root = arr[j]; // best root in loop
                for (int m = j; m<i; m++)
                {
                    if (fabs(kepler_equation(arr[m], e, M)) < fabs(kepler_equation(best_root, e, M)))
                        best_root =  arr[m];
                }

                return best_root;
            }
        }
    }

    return arr[MAX_ITER_NEWTON-1]; // in other cases return last approximation
}

double solve_kepler_eq(kepler_orbit* orbit, double M) {
    // Solve Kepler�s Equation for the eccentric anomaly
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

void kepler_to_cart(kepler_orbit* orbit, double grav_param,
                    double* pos, double* velo) {

    double M = orbit->M0;
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


void warp(double delta_t, double mass, kepler_orbit_denorm orbit_d, double* RA, double* dec) {
    double pos[3];
    double velo[3];
    kepler_orbit orbit;

    orbit_d.t0 = orbit_d.T0 + delta_t;
    normalize(&orbit_d, &orbit, R_BH_LY, mass);

    double grav = mass * G;
    double d = R_BH_LY;

    kepler_to_cart(&orbit, grav, pos, velo);

    pos[0] += 0;
    pos[1] += 0;
    pos[2] = d * LIGHT_YEAR;

    *RA = pos[1]/pos[2] * 180.0 * 3600.0 / M_PI;
    *dec = pos[0]/pos[2] * 180.0 * 3600.0 / M_PI;
}