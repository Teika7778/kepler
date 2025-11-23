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

#define EPS 1e6


void gauss_newton(double* parameters)
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

    // Перемеменная цикла Г-Н
    double* cur_val;

    // Начальное значение
    cur_val = parameters;

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
    // Стуктура для метода Рунге-Кутты 4
    rk4 rk_4 = {NULL, NULL, NULL, NULL, NULL};

    // Количество параметров
    int size = 7;

    // Массив векторов состояний численных производных
    double** deriv_state = (double**)malloc(sizeof(double*)*14);
    for (int j=0; j<14; j++)
        deriv_state[j] = (double*)malloc(sizeof(double)*STATE_SIZE_STAR);

    // Матрица AtWA
    double** AtWA = (double**)malloc(sizeof(double*)*size);
    for (int j=0; j<size; j++) AtWA[j] = (double*)malloc(sizeof(double)*size);

    // Вектор AtWr(betha)
    double AtWr[size];

    // Массив значений производных
    double deriv_val[14];

    while(++i != MAX_ITER_GAUSS_NEWTON)
    {

        double sum = 0;

        char buffer[256]; // Буфер для хранения строки

        // Заполнение нулями AtWA и AtWr(betha)
        for (int m=0; m<size; m++)
        {
            AtWr[m] = 0;
            for (int k=0; k<size; k++) AtWA[m][k] = 0;
        }

        // Цикл по звездвм
        for (size_t file_number=0; file_number<1; file_number++)
        {
            rewind(files[file_number]);

            // Начало интегрирования в первом наблюдении (ХАРДКОД)
            if (file_number == 0) previous_t = 2002.578;
            if (file_number == 1) previous_t = 2004.511;
            if (file_number == 2) previous_t = 2004.511;

            // Инициализация вектора системы текущими значениями параметров
            for(int j=0; j<STATE_SIZE_STAR; j++)
                x[j] = cur_val[j + STATE_SIZE_STAR*file_number];
            
            // Копируем вектор состояния в векотры производных
            for(int j=0; j<14; j++)
                memcpy(deriv_state[j], x, sizeof(double)*STATE_SIZE_STAR);

            // Добавляем и вычитаем eps
            for(int j=0; j<6; j++)
            {
                deriv_state[j*2][j] += std::abs(cur_val[j+STATE_SIZE_STAR*file_number] / EPS);
                deriv_state[j*2+1][j] -= std::abs(cur_val[j+STATE_SIZE_STAR*file_number] / EPS);
            }
                     
            while (fgets(buffer, sizeof(buffer), files[file_number]) != NULL)
            {

            // Чтение данных из файла
            sscanf(buffer, "%lf %lf %lf %lf %lf",
                   &t, &ra, &dec, &ra_err, &dec_err);

            // Численное интегирование:
            // Вектор системы
            wrap_integration(x, (t-previous_t)*365.*86400., cur_val[size-1], rk_4);

            // Интегрирование векторов производных не по массе
            for (int j=0; j<12; j++)
                wrap_integration(deriv_state[j], (t-previous_t)*365.*86400., cur_val[size-1], rk_4);

            // Интегрирование векторов производных по массе (требует eps в wrap_integration)
            wrap_integration(deriv_state[12], (t-previous_t)*365.*86400., cur_val[size-1]+(cur_val[size-1]/EPS), rk_4);
            wrap_integration(deriv_state[13], (t-previous_t)*365.*86400., cur_val[size-1]-(cur_val[size-1]/EPS), rk_4);

            // Вычисление невязки и проивзодных

            // Невязка
            g_i[0] = c/d* x[1]; //ra под 1
            g_i[1] = c/d* x[0];

            r_i[0] = g_i[0] - ra;
            r_i[1] = g_i[1] - dec;

            sum += pow(r_i[0], 2) / pow(ra_err, 2);
            sum += pow(r_i[1], 2) / pow(dec_err, 2);

            // Производные не по массе
            for(int j=0; j<6; j++)
            {
                // Центральные разности
                // По ra
                deriv_val[j*2] = 
                c/d*(deriv_state[j*2][1] - deriv_state[j*2+1][1])/(2*std::abs(cur_val[j+STATE_SIZE_STAR*file_number]/EPS));
                // по dec
                deriv_val[j*2+1]= 
                c/d*(deriv_state[j*2][0] - deriv_state[j*2+1][0])/(2*std::abs(cur_val[j+STATE_SIZE_STAR*file_number]/EPS));

                // Производные по z и v_z не нужно маштабировать
                if (j==2 or j==5)
                {
                    deriv_val[j*2] /= (c/d);
                    deriv_val[j*2+1] /= (c/d);
                }
            }

            //Производные по массе
            deriv_val[12] = c/d*(deriv_state[12][1] - deriv_state[13][1])/(2*std::abs(cur_val[size-1]/EPS));
            deriv_val[13] = c/d*(deriv_state[12][0] - deriv_state[13][0])/(2*std::abs(cur_val[size-1]/EPS));

            // Заполнение AtWr(betha)
            for(int j=0; j<size; j++)
            {
                AtWr[STATE_SIZE_STAR*file_number+j] += 
                (1.0/pow(ra_err, 2))*r_i[0]*deriv_val[j*2] + (1.0/pow(dec_err, 2))*r_i[1]*deriv_val[j*2+1];
            }

            // Заполнение AtWA

            // Матрица AtWA имеет блочный вид:

            // S_2    0    0     S_2M
            //  0   S_55   0     S_55M
            //  0    0    S_102  S_102M
            // S_2M S_55M S_102M S_2M+S_55M+S102_M

            for (int j=0; j<STATE_SIZE_STAR; j++)
            {
                // Заполнение блока конкретной звезды
                for (int k=0; k<6; k++)
                {
                    // Первая строка добавка по ra, вторая по dec
                    AtWA[j + STATE_SIZE_STAR*file_number][k + STATE_SIZE_STAR*file_number] += 
                    (1.0/pow(ra_err, 2)) * deriv_val[j*2] * deriv_val[k*2] +
                    (1.0/pow(dec_err, 2)) * deriv_val[j*2+1] * deriv_val[k*2+1];
                }

                // Заполнение правого столбца
                AtWA[j + STATE_SIZE_STAR*file_number][size-1] += 
                (1.0/pow(ra_err, 2)) * deriv_val[12] * deriv_val[j*2] +
                (1.0/pow(dec_err, 2)) * deriv_val[13] * deriv_val[j*2+1];
                
                // Заполнение нижней строки (Симметрия)
                AtWA[size-1][j + STATE_SIZE_STAR*file_number] = AtWA[j + STATE_SIZE_STAR*file_number][size-1];
                
                // Заполнение правого нижнего угла
                AtWA[size-1][size-1] += 
                1.0/pow(ra_err, 2) * deriv_val[12] * deriv_val[12] +
                1.0/pow(dec_err, 2) * deriv_val[13] * deriv_val[13];
            }

            previous_t = t;

            }
        }

        std::cout << "------------ITERATION " << i << " -----------------" << std::endl;
        std::cout << std::endl;

        std::cout << "ERROR SUM: " << sum << std::endl;
        std::cout << std::endl;

        for(int j=0; j<size; j++)
            printf("%.2e ", cur_val[j]);
        std::cout << std::endl;

        std::cout << std::endl;
        std::cout << std::endl;

        for(int j=0; j<size; j++)
            printf("%.2e ", AtWr[j]);
        std::cout << std::endl;

        std::cout << std::endl;
        std::cout << std::endl;

        for(int j=0; j<size; j++)
        {
            for(int k=0; k< size; k++)
                printf("%.2e ", AtWA[j][k]);
            std::cout << std::endl;
        }

        std::cout << std::endl;
        std::cout << std::endl;

        double w[size];

        solve_eq(AtWA, AtWr, size, w);

        for(int j=0; j<size; j++) cur_val[j] = cur_val[j] - w[j];

    }

    rk4Free(&rk_4);
    // Тут куча освобождений памяти

    fclose(files[0]);
    fclose(files[1]);
    fclose(files[2]);

}