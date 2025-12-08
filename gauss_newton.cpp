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


void find_inverse(double** A, double** res, int size)
{
    for (int i=0; i<size; i++)
    {
        double answ[size];
        double b[size];
        for(int j=0; j<size; j++)
        {
            if (j == i) b[j] = 1;
            else b[j] =0;
        }
        solve_eq(A, b, size, answ);

        for (int j=0; j<size; j++)
        {
            res[j][i] = answ[j];
        }
    }
}

void init_star(double* x, int star_number)
{
    double EPS = 5;
    switch (star_number)
    {
    case 0:
        x[0] = -1.34542432318288867188e+13 - -1.34542432318288867188e+13/ EPS;
        x[1] = 2.74369376973100244141e+12 + 2.74369376973100244141e+12 / EPS;
        x[2] = 1.17902665512916289062e+13 - 1.17902665512916289062e+13 / EPS;
        x[3] = 9.63030905550202805898e+03 + 9.63030905550202805898e+03 / EPS;
        x[4] = 2.38743901750857585284e+04 - 2.38743901750857585284e+04 / EPS;
        x[5] = 5.660631180e+03  + 5.660631180e+03/ EPS;
        break;
    case 1:
        x[0] = 4.0707e+12 - 4.0707e+12/ EPS;
        x[1] = 3.11869e+13 + 3.11869e+13 / EPS;
        x[2] = 2.54138e+12 - 2.54138e+12 / EPS;
        x[3] = 18926.6 + 18926.6 / EPS;
        x[4] = -2617.61 + -2617.61 / EPS;
        x[5] = 4412.44 - 4412.44/ EPS;
        break;
    case 2:
        x[0] = 3.06247e+13 - 3.06247e+13/ EPS;
        x[1] = -3.22222e+12 -3.22222e+12 / EPS;
        x[2] = 1.69223e+13 - 1.69223e+13 / EPS;
        x[3] = 1642.9 + 1642.9 / EPS;
        x[4] =-16531.6  -16531.6  / EPS;
        x[5] = -7379.3  + 7379.3 / EPS;
        break;

    default:
        x[0] = 0;
        x[1] = 0;
        x[2] = 0;
        x[3] = 0;
        x[4] = 0;
        x[5] = 0;
        break;
    }
}


void gauss_newton(double* parameters, int star_number, int* conditions)
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

    // Определяем размер в зависимости от количества параметров
    int size = 0;
    for(int j=0; j<6; j++)
        if (conditions[j] == 1) size += 1;

    int full_size = star_number*size + 1; 

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
        for (size_t file_number=0; file_number<star_number; file_number++)
        {
            rewind(files[file_number]);

            // Начало интегрирования в первом наблюдении (ХАРДКОД)
            if (file_number == 0) previous_t = 2002.578;
            if (file_number == 1) previous_t = 2004.511;
            if (file_number == 2) previous_t = 2004.511;

            init_star(x, file_number);

            int tmp = 0;

            // Инициализация вектора системы текущими значениями параметров
            for(int j=0; j<6; j++)
            {
                if (conditions[j] == 1)
                {
                    x[j] = cur_val[tmp + size*file_number];
                    tmp++;
                }
                
            }    

            // Инициализация авто производных
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

            for (int j=0; j<6; j++)
            {
                deriv[j*2] = x[12 + 6*1 + j];
                deriv[j*2+1] = x[12 + 6*0 + j];

                // Маштабирем наблюдаемые величины
                if (j == 0 or j == 1)
                {
                    deriv[j*2] *= c/d;
                    deriv[j*2+1] *= c/d;
                }
            }

            deriv[10] = c/d * x[7];
            deriv[11] = c/d * x[6];

            // Вычисление невязки и проивзодных

            // Невязка
            r_i[0] = c/d* x[1] - ra; // ra под 1
            r_i[1] = c/d* x[0] - dec; // dec под 0

            // Взвешенная сумма квадратов невязок
            sum += pow(r_i[0], 2) / pow(ra_err, 2);
            sum += pow(r_i[1], 2) / pow(dec_err, 2);

            // Техническая переменная для заполнения AtWr
            // Работает как второй счетчик цикла, который срабатывает
            // Только на тех значениях где conditions[j] == 1
            int t0 = 0;

            // Заполнение AtWr(betha)
            for(int j=0; j<6; j++)
            {
                if (conditions[j] == 1)
                {
                    AtWr[size*file_number+t0] += 
                    (1.0/pow(ra_err, 2))*r_i[0]*deriv[j*2] + (1.0/pow(dec_err, 2))*r_i[1]*deriv[j*2+1];
                    t0++;
                }
                
            }

            AtWr[full_size-1] += (1.0/pow(ra_err, 2))*r_i[0]*deriv[10] + (1.0/pow(dec_err, 2))*r_i[1]*deriv[11];
            

            // Заполнение AtWA

            // Матрица AtWA имеет блочный вид:

            // S_2    0    0     S_2M
            //  0   S_55   0     S_55M
            //  0    0    S_102  S_102M
            // S_2M S_55M S_102M S_2M+S_55M+S102_M

            // Аналогичные прошлой тех перменные
            int t1 = 0, t2=0;

            for (int j=0; j<6; j++)
            {
                if (conditions[j] == 1)
                {
                    t2 = 0;
                    // Заполнение блока конкретной звезды
                    for (int k=0; k<6; k++)
                    {
                        if (conditions[k] == 1)
                        {
                            // Первая строка добавка по ra, вторая по dec
                            AtWA[t1 + size*file_number][t2 + size*file_number] += 
                            (1.0/pow(ra_err, 2)) * deriv[j*2] * deriv[k*2] +
                            (1.0/pow(dec_err, 2)) * deriv[j*2+1] * deriv[k*2+1];
                            t2++;
                        }
                    }

                    // Заполнение правого столбца
                    AtWA[t1 + size*file_number][full_size-1] += 
                    (1.0/pow(ra_err, 2)) * deriv[10] * deriv[j*2] +
                    (1.0/pow(dec_err, 2)) * deriv[11] * deriv[j*2+1];

                    // Заполнение нижней строки (Симметрия)
                    AtWA[full_size-1][t1 + size*file_number] = AtWA[t1 + size*file_number][full_size-1];

                    t1++;
                }
            }

            // Заполнение правого нижнего угла
            AtWA[full_size-1][full_size-1] += 
            1.0/pow(ra_err, 2) * deriv[10] * deriv[10] +
            1.0/pow(dec_err, 2) * deriv[11] * deriv[11];


            previous_t = t;

            }
        }

        std::cout << "------------ITERATION " << i << " -----------------" << std::endl;
        std::cout << std::endl;

        std::cout << "ERROR SUM: " << sum << std::endl;
        std::cout << std::endl;

        for(int j=0; j<full_size; j++)
            printf("%.2e ", cur_val[j]);
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

        double w[full_size];

        solve_eq(AtWA, AtWr, full_size, w);

        for(int j=0; j<full_size; j++) cur_val[j] = cur_val[j] - w[j];

    }

    // Освобождение памяти
    rk4Free(&rk_4);

    for (int j=0; j<size; j++)
        free(AtWA[j]);
    free(AtWA);

    fclose(files[0]);
    fclose(files[1]);
    fclose(files[2]);

}