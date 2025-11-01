#include <cmath>
#include <iostream>

#include "struct.hpp"
#include "constans.hpp"
#include "diff.hpp"
#include "transform.hpp"

#include "gauss_newton.hpp"


void eval(kepler_orbit_denorm* stars, double M_bh)
{
    int i = 0;

    FILE* files[3];
    files[0] = fopen("data/s2.txt", "r");
    files[1] = fopen("data/s38.txt", "r");
    files[2] = fopen("data/s55.txt", "r");

    double numerator, denominator;
    double t, ra, dec, ra_err, dec_err;
    double g_i[2];
    double r_i[2];

    double d = R_BH_LY;

    numerator = 0;
    denominator = 0;
    double sum = 0;
    char buffer[256]; // Буфер для хранения строки
    for (size_t file_number=0; file_number<1; file_number++)
    {
        rewind(files[file_number]);
        while (fgets(buffer, sizeof(buffer), files[file_number]) != NULL) 
        {
        sscanf(buffer, "%lf %lf %lf %lf %lf", 
               &t, &ra, &dec, &ra_err, &dec_err);
        // count g_i (ra and dec from model)
        warp(t-stars[file_number].T0, M_bh, stars[file_number], &g_i[0], &g_i[1]);
        // count r_i
        r_i[0] = g_i[0] - ra;
        r_i[1] = g_i[1] - dec;
        sum += pow(r_i[0], 2) / pow(ra_err, 2);
        sum += pow(r_i[1], 2) / pow(dec_err, 2);
        }

        printf("MASS: %.20e, SS: %.20e\n", M_bh, sum);

    }

    fclose(files[0]);
    fclose(files[1]);
    fclose(files[2]);

}

double gauss_newton(kepler_orbit_denorm* stars, double M_bh)
{
    int i = 0;

    FILE* files[3];
    files[0] = fopen("data/s2.txt", "r");
    files[1] = fopen("data/s38.txt", "r");
    files[2] = fopen("data/s55.txt", "r");

    double arr[MAX_ITER_GAUSS_NEWTON];

    arr[0] = M_bh;

    double numerator, denominator;
    double t, ra, dec, ra_err, dec_err;
    double g_i[2];
    double r_i[2];
    double dr_i[2];

    double d = R_BH_LY;

    while(++i != MAX_ITER_GAUSS_NEWTON)
    {
        numerator = 0;
        denominator = 0;

        double sum = 0;

        char buffer[256]; // Буфер для хранения строки

        for (size_t file_number=0; file_number<3; file_number++)
        {
            rewind(files[file_number]);

            while (fgets(buffer, sizeof(buffer), files[file_number]) != NULL) 
            {

            sscanf(buffer, "%lf %lf %lf %lf %lf", 
                   &t, &ra, &dec, &ra_err, &dec_err);

            // count g_i (ra and dec from model)
            warp(t-stars[file_number].T0, arr[i-1], stars[file_number], &g_i[0], &g_i[1]);

            // count r_i

            r_i[0] = g_i[0] - ra;
            r_i[1] = g_i[1] - dec;

            sum += pow(r_i[0], 2) / pow(ra_err, 2);
            sum += pow(r_i[1], 2) / pow(dec_err, 2);

            //derivative_by_m(t, &stars[file_number], d*LIGHT_YEAR, arr[i-1], &dr_i[0], &dr_i[1]);

            count_diff(t, stars[file_number], arr[i-1], &dr_i[0], &dr_i[1], 1e24);

            numerator += ( 1./pow(ra_err, 2) ) * r_i[0] * dr_i[0];
            numerator += ( 1./pow(dec_err, 2) ) * r_i[1] * dr_i[1];

            denominator += ( 1./pow(ra_err, 2) ) * pow(dr_i[0], 2);
            denominator += ( 1./pow(dec_err, 2) ) * pow(dr_i[1], 2) ;

            }
        }

        printf("MASS: %.20e, SS: %.20e\n", arr[i-1], sum);

        arr[i] = arr[i-1] - numerator / denominator;

        if (std::isnan(arr[i]) || std::isinf(arr[i])) return arr[i];
    }

    fclose(files[0]);
    fclose(files[1]);
    fclose(files[2]);

    return arr[MAX_ITER_GAUSS_NEWTON-1];
}