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

// Количество параметров для одной звезды
#define SIZE 6
// Общее количество параметров
#define FULL_SIZE 19


void find_inverse(double** A, double** res, int size)
{
    for (int i=0; i<size; i++)
    {
        double answ[size];
        double b[size];
        for(int j=0; j<size; j++)
        {
            if (j == i) b[j] = 1;
            else b[j] = 0;
        }
        solve_eq(A, b, size, answ);

        for (int j=0; j<size; j++)
        {
            res[j][i] = answ[j];
        }
    }
}

void gauss_newton(double* parameters, int star_number, int* conditions, int GN_num_iter, double alpha)
{
    // Счетчик цикла Ньютона
    int i = 0;

    // Файлы данных наблюдений
    FILE* files[6];
    files[0] = fopen("data/s2.txt", "r");
    files[1] = fopen("data/s38.txt", "r");
    files[2] = fopen("data/s55.txt", "r");
    files[3] = fopen("data/s2_velocity.txt", "r");
    files[4] = fopen("data/s38_velocity.txt", "r");
    files[5] = fopen("data/s55_velocity.txt", "r");

    FILE* log = fopen("user_files/log.txt", "w");

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

    // Для скоростей
    double v_z, v_z_err, r_v_i;

    // Невязка
    double r_i[2];

    // Переменные для численного интегрирования:

    // Вектор состояния
    double x[STATE_SIZE_STAR_FULL];

    // Стуктура для метода Рунге-Кутты 4
    rk4 rk_4 = {NULL, NULL, NULL, NULL, NULL};

    // Значение малого возмущения
    double EPS = 1e8;

    // Матрица AtWA
    double** AtWA = (double**)malloc(sizeof(double*)*FULL_SIZE);
    for (int j=0; j<FULL_SIZE; j++) AtWA[j] = (double*)malloc(sizeof(double)*FULL_SIZE);

    // Вектор AtWr(betha)
    double AtWr[FULL_SIZE];

    // Предобуславливатель Якоби
    double Preconditioner[FULL_SIZE];

    // Копия кеплерова вектора системы для производных
    double* x_cep_copy = (double*) malloc(sizeof(double) * SIZE);

    // Массив производных
    double deriv[(SIZE+1) * 2];

    // Массив производных невязок vz по параметрам
    double deriv_vz[SIZE+1];

    while(++i != GN_num_iter)
    {

        double sum = 0;

        char buffer[256]; // Буфер для хранения строки

        // Заполнение нулями AtWA и AtWr(betha)
        for (int m=0; m<FULL_SIZE; m++)
        {
            AtWr[m] = 0;
            for (int k=0; k<FULL_SIZE; k++) AtWA[m][k] = 0;
        }

        // Цикл по звездвм
        for (size_t file_number=0; file_number<star_number+3; file_number++)
        {

            // Счетчик, дающий понять о какой звезде идет речь
            // при работе с ее радиальными скоростями
            int file_ra_dec = file_number % 3;

            rewind(files[file_number]);
            rewind(files[file_ra_dec]);

            // Начало интегрирования в первом наблюдении (ХАРДКОД)
            if (file_ra_dec == 0) previous_t = 2002.578;
            if (file_ra_dec == 1) previous_t = 2004.511;
            if (file_ra_dec == 2) previous_t = 2004.511;

            // Инициализация вектора системы текущими значениями параметров
            for(int j=0; j<6; j++)
                    x[j] = cur_val[j + SIZE*file_ra_dec];

            // Перевод из кеплеровых координат в декартовы
            star_state_to_dec(x, cur_val[FULL_SIZE-1], previous_t);

            // Инициализация авто производных (не используются)
            init_deriv(x);

            while (fgets(buffer, sizeof(buffer), files[file_number]) != NULL)
            {

            // Чтение данных из файла
            if (file_number > 2)
            {
                // Радиальные скорости
                sscanf(buffer, "%lf %lf %lf", &t, &v_z, &v_z_err);
            }
            else
            {
                // Ra dec
                sscanf(buffer, "%lf %lf %lf %lf %lf",
                   &t, &ra, &dec, &ra_err, &dec_err);
            }

            // Численное интегирование:
            // Вектор системы
            wrap_integration(x, (t-previous_t)*365.25*86400., cur_val[FULL_SIZE-1], rk_4);

            // Вычисление невязки и проивзодных

            // Невязка

            if (file_number > 2)
            {
                // Для радиальных скоростей
                r_i[0] = x[5]/1000 - v_z;
                r_i[1] = 0;
            }
            else{
                // Для ra и dec
                r_i[0] = c/d* x[1] - ra; // ra под 1
                r_i[1] = c/d* x[0] - dec; // dec под 0
            }

            if (file_number > 2)
            {
                sum += pow(r_i[0], 2) / pow(v_z_err, 2);
            }
            else{
                // Взвешенная сумма квадратов невязок
                sum += pow(r_i[0], 2) / pow(ra_err, 2);
                sum += pow(r_i[1], 2) / pow(dec_err, 2);
            }

            // Аналитические производные

            // Копируем текущие значения кеп. элемент
            for(int j=0; j<SIZE; j++)
                    x_cep_copy[j] = cur_val[j + SIZE*file_ra_dec];

            // Аналитические производные
            full_analytic(deriv, deriv_vz, x_cep_copy, cur_val[FULL_SIZE-1], t);
            

            // Маштабирование производных
            for(int j=0; j<SIZE; j++)
                    deriv_vz[j] /=  1000;
            deriv_vz[SIZE] /=  1000;

            // Отключение производных по параметрам, которые не определяются

            for (int j=0; j<SIZE; j++)
            {
                if (conditions[j] == 0)
                {
                    deriv[j*2] = 0;
                    deriv[j*2 + 1] = 0;
                    deriv_vz[j] = 0;
                }
            }


            // Заполнение AtWr(betha)
            for(int j=0; j<SIZE; j++)
            {
                if (file_number > 2)
                    AtWr[SIZE*file_ra_dec+j] +=
                (1.0/pow(v_z_err, 2))*r_i[0]*deriv_vz[j];
                else
                    AtWr[SIZE*file_ra_dec+j] +=
                (1.0/pow(ra_err, 2))*r_i[0]*deriv[j*2] + (1.0/pow(dec_err, 2))*r_i[1]*deriv[j*2+1];
            }

            if (file_number>2)
                AtWr[FULL_SIZE-1] += (1.0/pow(v_z_err, 2))*r_i[0]*deriv_vz[SIZE];
            else
                AtWr[FULL_SIZE-1] +=
            (1.0/pow(ra_err, 2))*r_i[0]*deriv[(SIZE+1)*2-2] + (1.0/pow(dec_err, 2))*r_i[1]*deriv[(SIZE+1)*2-1];


            // Заполнение AtWA

            // Матрица AtWA имеет блочный вид:

            // S_2    0    0     S_2M
            //  0   S_55   0     S_55M
            //  0    0    S_102  S_102M
            // S_2M S_55M S_102M S_2M+S_55M+S102_M

            for (int j=0; j<6; j++)
            {
                    // Заполнение блока конкретной звезды
                for (int k=0; k<6; k++)
                {
                    if (file_number > 2)
                        AtWA[j + SIZE*file_ra_dec][k + SIZE*file_ra_dec] +=
                    (1.0/pow(v_z_err, 2)) * deriv_vz[j] * deriv_vz[k];
                    else
                        AtWA[j + SIZE*file_ra_dec][k + SIZE*file_ra_dec] +=
                    (1.0/pow(ra_err, 2)) * deriv[j*2] * deriv[k*2] +
                    (1.0/pow(dec_err, 2)) * deriv[j*2+1] * deriv[k*2+1];
                }

                // Заполнение правого столбца
                if (file_number>2)
                    AtWA[j + SIZE*file_ra_dec][FULL_SIZE-1] +=
                (1.0/pow(v_z_err, 2)) * deriv_vz[SIZE] * deriv_vz[j];
                else
                    AtWA[j + SIZE*file_ra_dec][FULL_SIZE-1] +=
                (1.0/pow(ra_err, 2)) * deriv[(SIZE+1)*2-2] * deriv[j*2] +
                (1.0/pow(dec_err, 2)) * deriv[(SIZE+1)*2-1] * deriv[j*2+1];
                // Заполнение нижней строки (Симметрия)

                AtWA[FULL_SIZE-1][j + SIZE*file_ra_dec] = AtWA[j + SIZE*file_ra_dec][FULL_SIZE-1];
                
            }

            // Заполнение правого нижнего угла
            if (file_number>2)
                AtWA[FULL_SIZE-1][FULL_SIZE-1] +=
            1.0/pow(v_z_err, 2) * deriv_vz[SIZE] * deriv_vz[SIZE];
            else
                AtWA[FULL_SIZE-1][FULL_SIZE-1] +=
            1.0/pow(ra_err, 2) * deriv[(SIZE+1)*2-2] * deriv[(SIZE+1)*2-2] +
            1.0/pow(dec_err, 2) * deriv[(SIZE+1)*2-1] * deriv[(SIZE+1)*2-1];

            previous_t = t;

            }
        }

        // Устранение вырожденности для отключенных параметров
        for (int i = 0; i < star_number; i++) {
            for (int j = 0; j < 6; j++) {
                if (conditions[j] == 0) {
                    
                    AtWA[j + SIZE*i][j + SIZE*i] = 1.0; // Единица на диагонали
                    AtWr[j + SIZE*i] = 0.0;             // Ноль в векторе невязок
                }
            }
        }

        for(int i=0; i<star_number; i++)
        {
            for (int j=0; j<6; j++)
            {
                    Preconditioner[SIZE*i+j] = 1./sqrt(AtWA[j+ SIZE*i][j+ SIZE*i]);
            }
        }

        // Заполнение значений предобуславливателя (масса)

        Preconditioner[FULL_SIZE-1] = 1./sqrt(AtWA[FULL_SIZE-1][FULL_SIZE-1]);

        // Диагональное предобуславливание Якоби

        for(int i=0; i<FULL_SIZE; i++)
        {
            AtWr[i] *= Preconditioner[i];
            for(int j=0; j<FULL_SIZE; j++)
                AtWA[i][j] *= Preconditioner[i]*Preconditioner[j];
        }


        // --- ЗАПИСЬ В ЛОГ ФАЙЛ ---

        if (i == 1)
            std::cout << "Erorr sum on first iteration: " << sum << std::endl;
        if (i == GN_num_iter - 1)
            std::cout << "Erorr sum on iteration " << i << ": " << sum << std::endl;
        
        fprintf(log, "------------ITERATION %d -----------------\n\n", i);
        fprintf(log, "ERROR SUM: %g\n\n", sum);

        if (std::isnan(sum)) {
            // Если получили NaN, закрываем все файлы перед экстренным выходом
            fclose(log); 
            for (int f = 0; f < 6; f++) fclose(files[f]);
            return;
        }

        // Текущие значения параметров
        for(int j = 0; j < FULL_SIZE; j++) {
            fprintf(log, "%.2e ", cur_val[j]);
        }
        fprintf(log, "\n\n\n");

        // Вектор правой части (AtWr)
        for(int j = 0; j < FULL_SIZE; j++) {
            fprintf(log, "%.4e ", AtWr[j]);
        }
        fprintf(log, "\n\n\n");

        // Матрица системы (AtWA)
        for(int j = 0; j < FULL_SIZE; j++) {
            for(int k = 0; k < FULL_SIZE; k++) {
                fprintf(log, "%.2e ", AtWA[j][k]);
            }
            fprintf(log, "\n");
        }
        fprintf(log, "\n\n");
        
        // Чтобы данные сразу сохранялись на диск (полезно при падении программы)
        fflush(log); 
        
        // -------------------------

        double w[FULL_SIZE];
        solve_eq(AtWA, AtWr, FULL_SIZE, w);

        for(int j = 0; j < FULL_SIZE; j++) {
            cur_val[j] = cur_val[j] - alpha * Preconditioner[j] * w[j];
        }

    }

    // Освобождение памяти
    rk4Free(&rk_4);

    for (int j=0; j<SIZE; j++)
        free(AtWA[j]);
    free(AtWA);

    fclose(files[0]);
    fclose(files[1]);
    fclose(files[2]);
    fclose(files[3]);
    fclose(files[4]);
    fclose(files[5]);

    fclose(log);
}
