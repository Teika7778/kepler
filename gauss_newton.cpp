#include <cmath>
#include <iostream>
#include <string.h>

#include "struct.hpp"
#include "constans.hpp"
#include "diff.hpp"
#include "transform.hpp"
#include "integration.hpp"

#include "gauss_newton.hpp"

#define STEP_TO_PLOT_DIFF 10
#define STAR_TO_PLOT_DIFF 1 // s_38


double gauss_newton(kepler_orbit_denorm* stars, double M_bh, int steps, bool numerical)
{

    if (STEP_TO_PLOT_DIFF >= steps)
    {
        std::cout << "Слишком мало шагов для указанного STEP_TO_PLOT_DIF\n" << std::endl;
        return 0;
    }

    double d = (double) R_BH_LY * (double) LIGHT_YEAR;
    double c = 180 / M_PI * 3600;

    // Счетчик цикла Ньютона
    int i = 0;

    // Файлы данных наблюдений
    FILE* files[3];
    files[0] = fopen("data/s2.txt", "r");
    files[1] = fopen("data/s38.txt", "r");
    files[2] = fopen("data/s55.txt", "r");

    FILE* files_deriv[2];

    // Файл для записи результата
    if (numerical)
    {
        files_deriv[0] = fopen("plot_deriv/num_deriv_value.txt", "w+");
        files_deriv[1] = fopen("plot_deriv/num_deriv_error_sum.txt", "w+");
    }else{
        files_deriv[0] = fopen("plot_deriv/analytic_deriv_value.txt", "w+");
        files_deriv[1] = fopen("plot_deriv/analytic_deriv_error_sum.txt", "w+");
    }

    double current_value;

    // Начальное значение
    current_value = M_bh;

    // Числитель и знаменатель для вычисления нового шага:
    double numerator, denominator;

    // Переменные, считваемые из файла:
    double t, ra, dec, ra_err, dec_err;
    double previous_t;

    // Невязка
    double r_i[2];

    // Взвешенная сумма невязок:
    double sum;

    // Производная
    double dr_i[2];

    // Переменные для численного интегрирования:

    // Вектор состояния
    double x[STATE_SIZE_STAR];
    // Векторы состояния для орбит с инкрементом по массе
    double x_r[STATE_SIZE_STAR], x_l[STATE_SIZE_STAR];
    // Инкремент параметра (зависит от шага)
    double eps;
    // Стуктура для метода Рунге-Кутты 4
    rk4 rk_4 = {NULL, NULL, NULL, NULL, NULL};

    while(++i != steps)
    {
        numerator = 0, denominator = 0, sum = 0;

        eps = current_value / 1e8;

        char buffer[256]; // Буфер для хранения строки

        for (size_t file_number=0; file_number<3; file_number++)
        {
            rewind(files[file_number]);

            previous_t = stars[file_number].t0;

            init_star_state(x, stars[file_number], current_value); // init
            init_star_state(x_r, stars[file_number], current_value+eps); // init
            init_star_state(x_l, stars[file_number], current_value-eps); // init

            dr_i[0] = 0;
            dr_i[1] = 0;

            while (fgets(buffer, sizeof(buffer), files[file_number]) != NULL)
            {

            // Чтение данных из файла
            sscanf(buffer, "%lf %lf %lf %lf %lf",
                   &t, &ra, &dec, &ra_err, &dec_err);

            if (i == STEP_TO_PLOT_DIFF and file_number == STAR_TO_PLOT_DIFF)
                fprintf(files_deriv[0], "%.10e %.10e\n", dr_i[0], dr_i[1]);
            

            // Численное интегирование
            wrap_integration(x, (t-previous_t)*365.*86400., current_value, rk_4);
            wrap_integration(x_r, (t-previous_t)*365.*86400., current_value+eps, rk_4);
            wrap_integration(x_l, (t-previous_t)*365.*86400., current_value-eps, rk_4);

            if (numerical)
            {   // Численная
                dr_i[0] = c/d * (x_r[1] - x_l[1]) / (2*eps);
                dr_i[1] =  c/d * (x_r[0] - x_l[0]) / (2*eps);
            } else
            {   // Аналитическая
                dr_i[0] = c/d * x[7];
                dr_i[1] = c/d * x[6];
            }
            
            
            // Невязка
            r_i[0] = c/d* x[1] - ra;
            r_i[1] = c/d* x[0] - dec;

            // Взвешенная сумма невязок
            sum += pow(r_i[0], 2) / pow(ra_err, 2);
            sum += pow(r_i[1], 2) / pow(dec_err, 2);

            numerator += 1./pow(ra_err, 2)* r_i[0]*dr_i[0] +  1./pow(dec_err, 2)*r_i[1]*dr_i[1];
            denominator += 1./pow(ra_err, 2) * pow(dr_i[0], 2) + 1./pow(dec_err, 2)*pow(dr_i[1], 2);

            previous_t = t;

            }
        }

        fprintf(files_deriv[1], "%d %.10e\n", i, sum);

        current_value -=  numerator / denominator;

        if (std::isnan(current_value) || std::isinf(current_value))
        {
            fclose(files[0]);
            fclose(files[1]);
            fclose(files[2]);
            rk4Free(&rk_4);
            return current_value;
        }
    }

    rk4Free(&rk_4);

    fclose(files[0]);
    fclose(files[1]);
    fclose(files[2]);

    fclose(files_deriv[0]);
    fclose(files_deriv[1]);

    return current_value;
}
