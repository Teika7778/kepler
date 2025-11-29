#include <cmath>
#include <iostream>
#include <string.h>
#include <iomanip>

#include "struct.hpp"
#include "constans.hpp"
#include "diff.hpp"
#include "transform.hpp"
#include "integration.hpp"
#include "chol.hpp"
#include "helper_functions.hpp"

#include "gauss_newton.hpp"

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

    // Невязка
    double r_i[2];

    // Переменные для численного интегрирования

    // Вектор состояния
    double x[STATE_SIZE_STAR_FULL];
    // Стуктура для метода Рунге-Кутты 4
    rk4 rk_4 = {NULL, NULL, NULL, NULL, NULL};

    // Количество параметров сейчас ХАРДКОД для одной звезды
    int size = 6;
    int full_size = 16; 
    

    // Матрица AtWA
    double** AtWA = (double**)malloc(sizeof(double*)*full_size);
    for (int j=0; j<full_size; j++) AtWA[j] = (double*)malloc(sizeof(double)*full_size);

    // Вектор AtWr(betha)
    double AtWr[full_size];

    while(++i != MAX_ITER_GAUSS_NEWTON)
    {

        double sum = 0;

        char buffer[256]; // Буфер для хранения строки

        // Заполнение нулями AtWA и AtWr(betha)
        for (int m=0; m<full_size; m++)
        {
            AtWr[m] = 0;
            for (int k=0; k<full_size; k++) AtWA[m][k] = 0;
        }

        // Цикл по звездвм
        for (size_t file_number=0; file_number<3; file_number++)
        {
            rewind(files[file_number]);

            // Начало интегрирования в первом наблюдении (ХАРДКОД)
            if (file_number == 0) previous_t = 2002.578;
            if (file_number == 1) previous_t = 2004.511;
            if (file_number == 2) previous_t = 2004.511;

            // Инициализация вектора системы текущими значениями параметров
            for(int j=0; j<5; j++)
                x[j] = cur_val[j + 5*file_number];

            // Скорость по z не определяем
            x[5] = 5.660631180e+03;

            init_deriv(x);
                     
            while (fgets(buffer, sizeof(buffer), files[file_number]) != NULL)
            {

            // Чтение данных из файла
            sscanf(buffer, "%lf %lf %lf %lf %lf",
                   &t, &ra, &dec, &ra_err, &dec_err);

            // Численное интегирование:
            // Вектор системы
            wrap_integration(x, (t-previous_t)*365.*86400., cur_val[full_size-1], rk_4);

            // Массив производных
            double deriv[12];

            for (int j=0; j<size-1; j++)
            {
                deriv[j*2] = x[12 + 5*1 + j];
                deriv[j*2+1] = x[12 + 5*0 + j];

                // Маштабирем наблюдаемые величины
                if (j == 0 or j == 1)
                {
                    deriv[j*2] *= c/d;
                    deriv[j*2+1] *= c/d;
                }
            }

            deriv[10] = c/d * x[7];
            deriv[11] = c/d * x[6];

            //std::cout << "DERIVATIVES" << std::endl;
            //for (int j=0; j<6; j++)
            //    std::cout << deriv[2*j] <<  " " << deriv[2*j+1]  <<std::endl;
            //std::cout << std::endl;


            // Вычисление невязки и проивзодных

            // Невязка
            r_i[0] = c/d* x[1] - ra; // ra под 1
            r_i[1] = c/d* x[0] - dec; // dec под 0

            // Взвешенная сумма квадратов невязок
            sum += pow(r_i[0], 2) / pow(ra_err, 2);
            sum += pow(r_i[1], 2) / pow(dec_err, 2);

            // Заполнение AtWr(betha)
            for(int j=0; j<5; j++)
            {
                AtWr[5*file_number+j] += 
                (1.0/pow(ra_err, 2))*r_i[0]*deriv[j*2] + (1.0/pow(dec_err, 2))*r_i[1]*deriv[j*2+1];
            }
            

            // Заполнение AtWA

            // Матрица AtWA имеет блочный вид:

            // S_2    0    0     S_2M
            //  0   S_55   0     S_55M
            //  0    0    S_102  S_102M
            // S_2M S_55M S_102M S_2M+S_55M+S102_M

            for (int j=0; j<size-1; j++)
            {
                // Заполнение блока конкретной звезды
                for (int k=0; k<size-1; k++)
                {
                    // Первая строка добавка по ra, вторая по dec
                    AtWA[j + 5*file_number][k + 5*file_number] += 
                    (1.0/pow(ra_err, 2)) * deriv[j*2] * deriv[k*2] +
                    (1.0/pow(dec_err, 2)) * deriv[j*2+1] * deriv[k*2+1];
                }

                // Заполнение правого столбца
                AtWA[j + 5*file_number][full_size-1] += 
                (1.0/pow(ra_err, 2)) * deriv[10] * deriv[j*2] +
                (1.0/pow(dec_err, 2)) * deriv[11] * deriv[j*2+1];
                
                // Заполнение нижней строки (Симметрия)
                AtWA[full_size-1][j + 5*file_number] = AtWA[j + 5*file_number][full_size-1];
                
                // Заполнение правого нижнего угла
                AtWA[full_size-1][full_size-1] += 
                1.0/pow(ra_err, 2) * deriv[10] * deriv[10] +
                1.0/pow(dec_err, 2) * deriv[11] * deriv[11];
            }

            previous_t = t;

            }
        }

        std::cout << "------------ITERATION " << i << " -----------------" << std::endl;
        std::cout << std::endl;

        std::cout << "ERROR SUM: " << sum << std::endl;
        std::cout << std::endl;


        for(int j=0; j<full_size; j++)
            printf("%.6e ", cur_val[j]);
        std::cout << std::endl;

        

        std::cout << std::endl;
        std::cout << std::endl;

        for(int j=0; j<full_size; j++)
            printf("%.2e ", AtWr[j]);
        std::cout << std::endl;

        std::cout << std::endl;
        std::cout << std::endl;

        for(int j=0; j<full_size; j++)
        {
            for(int k=0; k< full_size; k++)
                printf("%.2e ", AtWA[j][k]);
            std::cout << std::endl;
        }

        std::cout << std::endl;
        std::cout << std::endl;

        /*
        

        double w[size];
        int num = 5;

        double* temp_matrix[num];
        for (int j=0; j<num; j++)
            temp_matrix[j] = (double*)malloc(sizeof(double)*num);

        double temp_vector[num];

        for (int j=0; j<num; j++)
        {
            for (int k=0; k<num; k++)
                temp_matrix[j][k] = AtWA[j][k];
            temp_vector[j] = AtWr[j];
        }

        

        solve_eq(temp_matrix, temp_vector, num, w);

        for (int j=0; j<num; j++)
            free(temp_matrix[j]);

        for(int j=num; j<7; j++)
            w[j] =0;
        w[5] = 0;

        for(int j=0; j<size; j++) cur_val[j] = cur_val[j] - w[j];
        */

        double w[full_size];

        solve_eq(AtWA, AtWr, full_size, w);

        for(int j=0; j<full_size; j++) cur_val[j] = cur_val[j] - w[j];

    }

    // Освобождение памяти
    rk4Free(&rk_4);

    //for (int j=0; j<14; j++)
    //    free(deriv_state[j]);
    for (int j=0; j<size; j++)
        free(AtWA[j]);
    free(AtWA);

    fclose(files[0]);
    fclose(files[1]);
    fclose(files[2]);

}