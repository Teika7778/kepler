#include <cmath>

#include "diff.hpp"
#include "constans.hpp"
#include "transform.hpp"

#define MAX_ITER_NEWTON 100

const double TWO_PI = 2.0 * M_PI;

double normalize_anomaly(double M) {
    double norm = fmod(M, TWO_PI);
    if (norm < 0) norm += TWO_PI;
    return norm;
}

double solve_kepler_eq(double e, double M) {
    // Solve Kepler�s Equation for the eccentric anomaly
    return newtons_method(M, e, M);
}

void derivative_by_m(double dt, kepler_orbit_denorm* denorm, double R_0, double m, double* ra, double* dec) {
    // -- CONST --
    double norm_a = (denorm->a * M_PI * (R_0 / 648000.0)) * LIGHT_YEAR;
    double norm_e = denorm->e;
    double norm_w = denorm->w * M_PI / 180.0;
    double norm_omega = denorm->omega * M_PI / 180.0;
    double norm_i = denorm->i * M_PI / 180.0;
    double delta_t = dt;

    // Матричные коэффициенты
    double A0 = cos(norm_w) * cos(norm_omega) - sin(norm_w) * cos(norm_i) * sin(norm_omega);
    double B0 = sin(norm_w) * cos(norm_omega) + cos(norm_w) * cos(norm_i) * sin(norm_omega);
    double A1 = cos(norm_w) * sin(norm_omega) + sin(norm_w) * cos(norm_i) * cos(norm_omega);
    double B1 = cos(norm_w) * cos(norm_i) * cos(norm_omega) - sin(norm_w) * sin(norm_omega);

    // Производная среднего движения по массе
    double u_by_m = 0.5 * sqrt(G / (m * pow(norm_a, 3)));

    double M_by_m = delta_t * u_by_m;

    double M = fmod(sqrt( (m/pow(norm_a, 3)) * G) *delta_t, 2*M_PI);
    double E = solve_kepler_eq(norm_e, M);

    double E_by_m = M_by_m / (1 - norm_e * cos(E));
    double D_by_m = norm_a * norm_e * sin(E) * E_by_m;
    double v_by_m = E_by_m * sqrt(1 - pow(norm_e, 2))/(1 - norm_e * cos(E));

    double v = 2*atan2(
                sqrt(1+norm_e) * sin(E/2),
                sqrt(1-norm_e) * cos(E/2)
            );
    double D = norm_a * (1 - norm_e * cos(E));

    double x_by_m = D_by_m * cos(v) - D * sin(v) * v_by_m;
    double y_by_m = D_by_m * sin(v) + D * cos(v) * v_by_m;

    double pos1_by_m = x_by_m * A0 - y_by_m * B0;
    double pos0_by_m = x_by_m * A1 - y_by_m * B1;

    double pos2 = (double) R_BH_LY * (double) LIGHT_YEAR;

    *ra = (180.0 * 3600.0) / (pos2 * M_PI) * pos1_by_m;
    *dec = (180.0 * 3600.0) / (pos2 * M_PI) * pos0_by_m;
}



void count_diff(double dt, kepler_orbit_denorm denorm, double m, double* ra, double* dec, double eps)
{

    double coords_1[2];
    double coords_2[2];

    warp(dt - denorm.T0, m+eps, denorm, &coords_1[0], &coords_1[1]);

    warp(dt - denorm.T0, m-eps, denorm, &coords_2[0], &coords_2[1]);

    *ra = ( coords_1[0] - coords_2[0] ) / (2*eps);

    *dec = ( coords_1[1] - coords_2[1] ) / (2*eps);

}

void full_analytic(double* deriv_vec, double* params, double t) {
    // printf("t = %.4e\n", t);
    t *= 365.25 * 86400.0;
    double d = (double) R_BH_LY * (double) LIGHT_YEAR;
    double c = 180.0 / M_PI * 3600.0;

    for (int j = 0; j < 14; j++) {
        deriv_vec[j] = 0;
    }

    /*
    for (int j = 0; j < 7; j++) {
        printf("%.4e \n", params[j]);
    }
    printf("\n");
    */

    double a = params[0];
    double e = params[1];
    double w = params[2];
    double O = params[3];
    double i = params[4];
    double Tp = params[5];
    double m = params[6];

    double grav = G * m; // grav param
    double dgrav_dm = G;

    double na = d * a / c; //a srcsec --> meters
    double dna_da = d / c;

    double Tps = 365.25 * 86400.0 * Tp; //Tp yr --> sec
    double dTps_dTp = 365.25 * 86400.0;

    // deg --> rad
    double dangle = M_PI / 180.0;
    double wrad = w * dangle;
    double Orad = O * dangle;
    double irad = i * dangle;

    // printf("%.4e %.4e %.4e %.4e %.4e %.4e %.4e \n\n", na, e, wrad, Orad, irad, Tps, grav);

    double n = sqrt(grav / pow(na, 3));
    double dn_da = -1.5 * sqrt(grav / pow(na,5)) * dna_da;
    double dn_dm = dgrav_dm / (2 * sqrt(grav * pow(na, 3)));

    double dt = t - Tps;
    double ddt_dTp = -1 * dTps_dTp;

    double M = normalize_anomaly(n*dt);
    double dM_da = dt * dn_da;
    double dM_dm = dt * dn_dm;
    double dM_dTp = n * ddt_dTp;

    double E = solve_kepler_eq(e, M);
    double dE_de = sin(E) / (1 - e*cos(E));
    double dE_da = dM_da / (1 - e*cos(E));
    double dE_dm = dM_dm / (1 - e*cos(E));
    double dE_dTp = dM_dTp / (1 - e*cos(E));

    double v = 2*atan2(
                sqrt(1+e) * sin(E/2),
                sqrt(1-e) * cos(E/2)
            );
    double dv_de = ( (pow(e,2) - 1)*dE_de - sin(E) )  /  ( sqrt(1-pow(e,2))*(e*cos(E) - 1) );
    double dv_da = sqrt(1-pow(e,2)) * dE_da / (1 - e*cos(E));
    double dv_dm = sqrt(1-pow(e,2)) * dE_dm / (1 - e*cos(E));
    double dv_dTp = sqrt(1-pow(e,2)) * dE_dTp / (1 - e*cos(E));

    double rc = na*(1-e*cos(E));
    double drc_de = na*e*dE_de*sin(E) - na*cos(E);
    double drc_da = na*e*dE_da*sin(E) + dna_da * (1 - e*cos(E));
    double drc_dm = e * na * sin(E) * dE_dm;
    double drc_dTp = e * na * sin(E) * dE_dTp;

    double rx = rc * cos(v);
    double ry = rc * sin(v);
    double drx_de = drc_de * cos(v) - rc * sin(v) * dv_de;
    double dry_de = drc_de * sin(v) + rc * cos(v) * dv_de;
    double drx_da = drc_da * cos(v) - rc * sin(v) * dv_da;
    double dry_da = drc_da * sin(v) + rc * cos(v) * dv_da;
    double drx_dm = drc_dm * cos(v) - rc * sin(v) * dv_dm;
    double dry_dm = drc_dm * sin(v) + rc * cos(v) * dv_dm;
    double drx_dTp = drc_dTp * cos(v) - rc * sin(v) * dv_dTp;
    double dry_dTp = drc_dTp * sin(v) + rc * cos(v) * dv_dTp;


    double rrx = rx * (cos(wrad) * cos(Orad) - sin(wrad) * cos(irad) * sin(Orad))
               - ry * (sin(wrad) * cos(Orad) + cos(wrad) * cos(irad) * sin(Orad));

    double rry = rx * (cos(wrad) * sin(Orad) + sin(wrad) * cos(irad) * cos(Orad))
               + ry * (cos(wrad) * cos(irad) * cos(Orad) - sin(wrad) * sin(Orad));

    double drrx_de = drx_de * (cos(wrad) * cos(Orad) - sin(wrad) * cos(irad) * sin(Orad))
                  - dry_de * (sin(wrad) * cos(Orad) + cos(wrad) * cos(irad) * sin(Orad));

    double drry_de = drx_de * (cos(wrad) * sin(Orad) + sin(wrad) * cos(irad) * cos(Orad))
                  + dry_de * (cos(wrad) * cos(irad) * cos(Orad) - sin(wrad) * sin(Orad));

    double drrx_da = drx_da * (cos(wrad) * cos(Orad) - sin(wrad) * cos(irad) * sin(Orad))
                  - dry_da * (sin(wrad) * cos(Orad) + cos(wrad) * cos(irad) * sin(Orad));

    double drry_da = drx_da * (cos(wrad) * sin(Orad) + sin(wrad) * cos(irad) * cos(Orad))
                  + dry_da * (cos(wrad) * cos(irad) * cos(Orad) - sin(wrad) * sin(Orad));

    double drrx_dm = drx_dm * (cos(wrad) * cos(Orad) - sin(wrad) * cos(irad) * sin(Orad))
                  - dry_dm * (sin(wrad) * cos(Orad) + cos(wrad) * cos(irad) * sin(Orad));

    double drry_dm = drx_dm * (cos(wrad) * sin(Orad) + sin(wrad) * cos(irad) * cos(Orad))
                  + dry_dm * (cos(wrad) * cos(irad) * cos(Orad) - sin(wrad) * sin(Orad));

    double drrx_dTp = drx_dTp * (cos(wrad) * cos(Orad) - sin(wrad) * cos(irad) * sin(Orad))
                   - dry_dTp * (sin(wrad) * cos(Orad) + cos(wrad) * cos(irad) * sin(Orad));

    double drry_dTp = drx_dTp * (cos(wrad) * sin(Orad) + sin(wrad) * cos(irad) * cos(Orad))
                   + dry_dTp * (cos(wrad) * cos(irad) * cos(Orad) - sin(wrad) * sin(Orad));

    /*
    double rrx = rx * (cos(wrad) * cos(Orad) - sin(wrad) * cos(irad) * sin(Orad))
               - ry * (sin(wrad) * cos(Orad) + cos(wrad) * cos(irad) * sin(Orad));

    double rry = rx * (cos(wrad) * sin(Orad) + sin(wrad) * cos(irad) * cos(Orad))
               + ry * (cos(wrad) * cos(irad) * cos(Orad) - sin(wrad) * sin(Orad));
    */

    // ПРОИЗВОДНЫЕ ДОПИСАТЬ
    double drrx_dw = rx * (-1 * sin(wrad) * dangle * cos(Orad) - cos(wrad) * dangle * cos(irad) * sin(Orad))
                   - ry * (cos(wrad) * dangle * cos(Orad) - sin(wrad) * dangle * cos(irad) * sin(Orad));

    double drry_dw = rx * (-1 * sin(wrad) * dangle * sin(Orad) + cos(wrad) * dangle * cos(irad) * cos(Orad))
                   + ry * (-1 * sin(wrad) * dangle * cos(irad) * cos(Orad) - cos(wrad) * dangle * sin(Orad));

    double drrx_dO = rx * (-1 * cos(wrad) * sin(Orad) * dangle - sin(wrad) * cos(irad) * cos(Orad) * dangle)
                   - ry * (-1 * sin(wrad) * sin(Orad) * dangle + cos(wrad) * cos(irad) * cos(Orad) * dangle);

    double drry_dO = rx * (cos(wrad) * cos(Orad) * dangle - sin(wrad) * cos(irad) * sin(Orad) * dangle)
                   + ry * (-1 * cos(wrad) * cos(irad) * sin(Orad) * dangle - sin(wrad) * cos(Orad) * dangle);

    double drrx_di = rx * (sin(wrad) * sin(irad) * dangle * sin(Orad))
                   - ry * (-1 * cos(wrad) * sin(irad) * dangle * sin(Orad));

    double drry_di = rx * (-1 * sin(wrad) * sin(irad) * dangle * cos(Orad))
                   + ry * (-1 * cos(wrad) * sin(irad) * dangle * cos(Orad));

    double DEC = c/d * rrx;
    double RA = c/d * rry;

    double dDEC_de = c/d * drrx_de;
    double dDEC_da = c/d * drrx_da;
    double dDEC_dm = c/d * drrx_dm;
    double dDEC_dw = c/d * drrx_dw;
    double dDEC_dO = c/d * drrx_dO;
    double dDEC_di = c/d * drrx_di;
    double dDEC_dTp = c/d * drrx_dTp;

    double dRA_de = c/d * drry_de;
    double dRA_da = c/d * drry_da;
    double dRA_dm = c/d * drry_dm;
    double dRA_dw = c/d * drry_dw;
    double dRA_dO = c/d * drry_dO;
    double dRA_di = c/d * drry_di;
    double dRA_dTp = c/d * drry_dTp;
    /*
    three_stars[0] = 0.126;  // a
    three_stars[1] = 0.884;  // e
    three_stars[2] = 71.36;   // w
    three_stars[3] = 234.50;  // omega
    three_stars[4] = 136.78;  // i
    three_stars[5] = 2002.32; // T0
    */
    deriv_vec[0] = dRA_da;
    deriv_vec[1] = dDEC_da;
    deriv_vec[2] = dRA_de;
    deriv_vec[3] = dDEC_de;
    deriv_vec[4] = dRA_dw;
    deriv_vec[5] = dDEC_dw;
    deriv_vec[6] = dRA_dO;
    deriv_vec[7] = dDEC_dO;
    deriv_vec[8] = dRA_di;
    deriv_vec[9] = dDEC_di;
    deriv_vec[10] = dRA_dTp;
    deriv_vec[11] = dDEC_dTp;
    deriv_vec[12] = dRA_dm;
    deriv_vec[13] = dDEC_dm;
}
