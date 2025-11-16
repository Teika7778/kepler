#include <cmath>
#include <iostream>
#include <string.h>

#include "struct.hpp"
#include "constans.hpp"
#include "diff.hpp"
#include "transform.hpp"
#include "integration.hpp"

#include "gauss_newton.hpp"


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


double gauss_newton_2(kepler_orbit_denorm* stars, double M_bh)
{
    // Счетчик цикла Ньютона
    int i = 0;

    // Файлы данных наблюдений
    FILE* files[3];
    files[0] = fopen("data/s2.txt", "r");
    files[1] = fopen("data/s38.txt", "r");
    files[2] = fopen("data/s55.txt", "r");

    // Массив значений для метода Ньютона
    double arr[MAX_ITER_GAUSS_NEWTON];

    // Начальное значение
    arr[0] = M_bh;

    // Числитель и знаменатель для вычисления нового шага
    double numerator, denominator;
    // Переменные, считваемые из файла
    double t, ra, dec, ra_err, dec_err;
    double previous_t;
    // Рассчитанные значения на i-ом измерении
    double g_i[2];
    // Неувязка
    double r_i[2];
    // Рассчитанная производная
    double dr_i[2];

    // Переменные для численного интегрирования

    // Вектор состояния
    double x[STATE_SIZE_STAR];
    // Вектор состояния производных
    double deriv[STATE_SIZE_DERIV];
    // Стуктура для метода Рунге-Кутты 4
    rk4 rk_4 = {NULL, NULL, NULL, NULL, NULL};
    // Занчения координат (для передачи в производную)
    double x_for_deriv[1], y_for_deriv[1], z_for_deriv[1];
    double* array_for_deriv[] = {x_for_deriv, y_for_deriv, z_for_deriv};

    // Переменные для метода центральных разностей


    while(++i != MAX_ITER_GAUSS_NEWTON)
    {
        numerator = 0;
        denominator = 0;

        double sum = 0;

        char buffer[256]; // Буфер для хранения строки

        for (size_t file_number=0; file_number<2; file_number++)
        {
            rewind(files[file_number]);

            previous_t = stars[file_number].t0;

            init_star_state(x, stars[file_number], arr[i-1]); // init

            for (int j = 0; j < STATE_SIZE_DERIV; j++){
                deriv[j] = 0;
            }

            while (fgets(buffer, sizeof(buffer), files[file_number]) != NULL)
            {

            // Чтение данных из файла
            sscanf(buffer, "%lf %lf %lf %lf %lf",
                   &t, &ra, &dec, &ra_err, &dec_err);

            // Численное интегирование

            if (t < previous_t)
            {
                wrap_integration(x, deriv, (t-previous_t)*365.*86400., arr[i-1], rk_4, array_for_deriv);
            } else
            {
                wrap_integration(x, deriv, (t-previous_t)*365.*86400., arr[i-1], rk_4, array_for_deriv);
            }


            g_i[0] = x[1]; //ra под 1
            g_i[1] = x[0];


            double d = (double) R_BH_LY * (double) LIGHT_YEAR;
            double c = 180 / M_PI * 3600;

            dr_i[0] = deriv[1] * d / c;
            dr_i[1] = deriv[0] * d / c;
            // count r_i


            r_i[0] = g_i[0] - ra;
            r_i[1] = g_i[1] - dec;

            sum += pow(r_i[0], 2) / pow(ra_err, 2);
            sum += pow(r_i[1], 2) / pow(dec_err, 2);

            numerator += ( 1./pow(ra_err, 2) ) * r_i[0] * dr_i[0];
            numerator += ( 1./pow(dec_err, 2) ) * r_i[1] * dr_i[1];

            denominator += ( 1./pow(ra_err, 2) ) * pow(dr_i[0], 2);
            denominator += ( 1./pow(dec_err, 2) ) * pow(dr_i[1], 2);

            previous_t = t;

            }
        }

        printf("MASS: %.20e, SS: %.20e\n", arr[i-1], sum);

        arr[i] = arr[i-1] - numerator / denominator;

        if (std::isnan(arr[i]) || std::isinf(arr[i]))
        {
            fclose(files[0]);
            fclose(files[1]);
            fclose(files[2]);
            rk4Free(&rk_4);
            printf("MASS: %.20e, SS: %.20e\n", arr[i], sum);
            return arr[i];
        }

    }

    rk4Free(&rk_4);

    fclose(files[0]);
    fclose(files[1]);
    fclose(files[2]);

    return arr[MAX_ITER_GAUSS_NEWTON-1];
}
