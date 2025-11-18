#include <cmath>
#include <iostream>
#include <string.h>

#include "struct.hpp"
#include "constans.hpp"
#include "diff.hpp"
#include "transform.hpp"
#include "integration.hpp"
#include "chol.hpp"
#include "helper_functions.hpp"

#include "gauss_newton.hpp"


double gauss_newton(kepler_orbit_denorm* stars, double M_bh)
{

    double d = (double) R_BH_LY * (double) LIGHT_YEAR;
    double c = 180 / M_PI * 3600;

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
    // Невязка
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

        for (size_t file_number=0; file_number<3; file_number++)
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
            wrap_integration(x, deriv, (t-previous_t)*365.*86400., arr[i-1], rk_4, array_for_deriv);

            g_i[0] = x[1] * c / d; //ra под 1
            g_i[1] = x[0] * c / d;

            dr_i[0] = deriv[1] * c / d;
            dr_i[1] = deriv[0] * c / d;
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

        if (i == MAX_ITER_GAUSS_NEWTON-1)
            printf("MASS: %.20e, SS: %.20e\n", arr[i], sum);
    }

    

    rk4Free(&rk_4);

    fclose(files[0]);
    fclose(files[1]);
    fclose(files[2]);

    return arr[MAX_ITER_GAUSS_NEWTON-1];
}


double gauss_newton_3(double* parameters)
{
    // Счетчик цикла Ньютона
    int i = 0;

    // Файлы данных наблюдений
    FILE* files[3];
    files[0] = fopen("data/s2.txt", "r");
    files[1] = fopen("data/s38.txt", "r");
    files[2] = fopen("data/s55.txt", "r");

    // Переменные для перевода метров в ra. и dec.
    double d = (double) R_BH_LY * (double) LIGHT_YEAR;
    double c = 180 / M_PI * 3600;

    // Массив значений для метода Ньютона
    double* arr[MAX_ITER_GAUSS_NEWTON];

    // Начальное значение
    arr[0] = parameters;

    // Переменные, считваемые из файла
    double t, ra, dec, ra_err, dec_err;
    double previous_t;

    // Рассчитанные значения на i-ом измерении
    double g_i[2];
    // Невязка
    double r_i[2];

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

    // Количество параметров
    int size = sizeof(parameters) / sizeof(double);

    // Матрица AtWA
    double** AtWA = (double**)malloc(sizeof(double*)*size);
    for (int j=0; j<size; j++) AtWA[j] = (double*)malloc(sizeof(double)*size);

    // Вектор AtWr(betha)
    double* AtWr = (double*)malloc(sizeof(double)*size);

    // Численные производные
    double x_r_m[STATE_SIZE_STAR], x_l_m[STATE_SIZE_STAR]; // По массе
    double x_r_x0[STATE_SIZE_STAR], x_l_x0[STATE_SIZE_STAR]; // По x_0
    double x_r_y0[STATE_SIZE_STAR], x_l_y0[STATE_SIZE_STAR]; // По y_0
    double x_r_z0[STATE_SIZE_STAR], x_l_z0[STATE_SIZE_STAR]; // По z_0
    double x_r_v_x0[STATE_SIZE_STAR], x_l_v_x0[STATE_SIZE_STAR]; // По v_x0
    double x_r_v_y0[STATE_SIZE_STAR], x_l_v_y0[STATE_SIZE_STAR]; // По v_y0
    double x_r_v_z0[STATE_SIZE_STAR], x_l_v_z0[STATE_SIZE_STAR]; // По v_z0

    // Подходящие epsilon
    double e_m = 1e30;
    double e_x0 = 1e10, e_y0 = 1e10, e_z0 = 1e10;
    double e_v_x0 = 1e10, e_v_y0 = 1e10, e_v_z0 = 1e10;

    // Массивы для производных

    double dx_dm[2], dx_dx0[2], dx_dy0[2], dx_dz0[2], dx_dv_x0[2], dx_dv_y0[2], dx_dv_z0[2]; 

    while(++i != MAX_ITER_GAUSS_NEWTON)
    {

        char buffer[256]; // Буфер для хранения строки

        // Заполнение нулями
        for (int m=0; m<size; m++)
        {
            AtWr[m] = 0;
            for (int k=0; k<size; k++) AtWA[m][k] = 0;
        }

        for (size_t file_number=0; file_number<3; file_number++)
        {
            rewind(files[file_number]);

            // Начало интегрирования в первом наблюдении (ХАРДКОД)
            if (file_number == 0) previous_t = 2002.578;
            if (file_number == 1) previous_t = 2004.511;
            if (file_number == 2) previous_t = 2004.511;

            // Инициализация текущими значениями параметров
            for(int j=0; j<STATE_SIZE_STAR; j++)
            {
                // вектор системы
                x[j] = arr[i-1][j + STATE_SIZE_STAR*file_number];
                // векотры производных
                // масса 
                x_r_m[j] = arr[i-1][j + STATE_SIZE_STAR*file_number];
                x_l_m[j] = arr[i-1][j + STATE_SIZE_STAR*file_number];
                // x_0
                x_r_x0[j] = arr[i-1][j + STATE_SIZE_STAR*file_number];
                x_l_x0[j] = arr[i-1][j + STATE_SIZE_STAR*file_number];
                // y_0
                x_r_y0[j] = arr[i-1][j + STATE_SIZE_STAR*file_number];
                x_l_y0[j] = arr[i-1][j + STATE_SIZE_STAR*file_number];
                // z_0
                x_r_z0[j] = arr[i-1][j + STATE_SIZE_STAR*file_number];
                x_l_z0[j] = arr[i-1][j + STATE_SIZE_STAR*file_number];
                // v_x0
                x_r_v_x0[j] = arr[i-1][j + STATE_SIZE_STAR*file_number];
                x_l_v_x0[j] = arr[i-1][j + STATE_SIZE_STAR*file_number];
                // v_y0
                x_r_v_y0[j] = arr[i-1][j + STATE_SIZE_STAR*file_number];
                x_l_v_y0[j] = arr[i-1][j + STATE_SIZE_STAR*file_number];
                // v_z0
                x_r_v_z0[j] = arr[i-1][j + STATE_SIZE_STAR*file_number];
                x_l_v_z0[j] = arr[i-1][j + STATE_SIZE_STAR*file_number];
            }

            // Добавление eps
            x_r_x0[0] += e_x0;
            x_l_x0[0] -= e_x0;

            x_r_y0[0] += e_y0;
            x_l_y0[0] -= e_y0;

            x_r_z0[0] += e_z0;
            x_l_z0[0] -= e_z0;

            x_r_v_x0[0] += e_v_x0;
            x_l_v_x0[0] -= e_v_x0;

            x_r_v_y0[0] += e_v_y0;
            x_l_v_y0[0] -= e_v_y0;

            x_r_v_z0[0] += e_v_z0;
            x_l_v_z0[0] -= e_v_z0;
                     

            while (fgets(buffer, sizeof(buffer), files[file_number]) != NULL)
            {

            // Чтение данных из файла
            sscanf(buffer, "%lf %lf %lf %lf %lf",
                   &t, &ra, &dec, &ra_err, &dec_err);

            // Численное интегирование
            // Вектор системы
            wrap_integration(x, deriv, (t-previous_t)*365.*86400., arr[i-1][size-1], rk_4, array_for_deriv);

            // Численные производные
            wrap_integration(x_r_m, deriv, (t-previous_t)*365.*86400., arr[i-1][size-1]+e_m, rk_4, array_for_deriv);
            wrap_integration(x_l_m, deriv, (t-previous_t)*365.*86400., arr[i-1][size-1]-e_m, rk_4, array_for_deriv);

            wrap_integration(x_r_x0, deriv, (t-previous_t)*365.*86400., arr[i-1][size-1], rk_4, array_for_deriv);
            wrap_integration(x_l_x0, deriv, (t-previous_t)*365.*86400., arr[i-1][size-1], rk_4, array_for_deriv);

            wrap_integration(x_r_y0, deriv, (t-previous_t)*365.*86400., arr[i-1][size-1], rk_4, array_for_deriv);
            wrap_integration(x_l_y0, deriv, (t-previous_t)*365.*86400., arr[i-1][size-1], rk_4, array_for_deriv);

            wrap_integration(x_r_z0, deriv, (t-previous_t)*365.*86400., arr[i-1][size-1], rk_4, array_for_deriv);
            wrap_integration(x_l_z0, deriv, (t-previous_t)*365.*86400., arr[i-1][size-1], rk_4, array_for_deriv);

            wrap_integration(x_r_v_x0, deriv, (t-previous_t)*365.*86400., arr[i-1][size-1], rk_4, array_for_deriv);
            wrap_integration(x_l_v_x0, deriv, (t-previous_t)*365.*86400., arr[i-1][size-1], rk_4, array_for_deriv);

            wrap_integration(x_r_v_z0, deriv, (t-previous_t)*365.*86400., arr[i-1][size-1], rk_4, array_for_deriv);
            wrap_integration(x_l_v_z0, deriv, (t-previous_t)*365.*86400., arr[i-1][size-1], rk_4, array_for_deriv);

            // Вычисление невязки и проивзодных

            // Невязка
            g_i[0] = x[1] * d/c; //ra под 1
            g_i[1] = x[0] * d/c;

            r_i[0] = g_i[0] - ra;
            r_i[1] = g_i[1] - dec;

            // Производные

            dx_dm[0] = (x_r_m[1] - x_l_m[1]) / 2*e_m;

            dr_i[0] = deriv[1] * d / c;
            dr_i[1] = deriv[0] * d / c;

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