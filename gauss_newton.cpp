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
            else b[j] = 0;
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
    switch (star_number)
    {
    case 0:
        x[0] = 0.126;  // a
        x[1] = 0.884;  // e
        x[2] = 71.36;   // w
        x[3] = 234.50;  // omega
        x[4] = 136.78;  // i
        x[5] = 2002.32; // T0
        break;
    case 1:
        x[0] = 0.140;  // a
        x[1] = 0.818;  // e
        x[2] = 18.4;   // w
        x[3] = 101.8;  // omega
        x[4] = 166.22;  // i
        x[5] = 2003.30; // T0
        break;
    case 2:
        x[0] = 0.109;  // a
        x[1] = 0.74;  // e
        x[2] = 133.5;   // w
        x[3] = 129.9;  // omega
        x[4] = 141.7;  // i
        x[5] = 2009.31; // T0
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
    // Все возможнные производные (потом будут перенесены только нужные)
    double* deriv_a = (double*) malloc(sizeof(double) * 14);
    double* deriv_in = (double*) malloc(sizeof(double) * 7);
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

    // кол-во орбит для центральных разностей
    int deriv_arr_size = (size + 1) * 2;

    // Значение малого возмущения
    double EPS = 1e2;

    // Массив векторов состояний численных производных
    double* deriv_state[deriv_arr_size];
    for (int j=0; j<deriv_arr_size; j++)
        deriv_state[j] = (double*)malloc(sizeof(double)*STATE_SIZE_STAR_FULL);

    // Матрица AtWA
    double** AtWA = (double**)malloc(sizeof(double*)*full_size);
    for (int j=0; j<full_size; j++) AtWA[j] = (double*)malloc(sizeof(double)*full_size);

    // Вектор AtWr(betha)
    double AtWr[full_size];

    // Предобуславливатель Якоби
    double Preconditioner[full_size];

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

            // Копируем вектор состояния в векотры производных
            for(int j=0; j<deriv_arr_size; j++)
                memcpy(deriv_state[j], x, sizeof(double)*STATE_SIZE_STAR_FULL);

            tmp = 0;
            // Добавляем и вычитаем eps (кроме производных по массе)
            for(int j=0; j<6; j++)
            {
                if (conditions[j] == 1)
                {
                    deriv_state[tmp*2][j] += std::abs(cur_val[tmp+size*file_number] / EPS);
                    deriv_state[tmp*2+1][j] -= std::abs(cur_val[tmp+size*file_number] / EPS);
                    tmp++;
                }
            }

            //std::cout << "Первая T0: " << deriv_state[0][5] << "\n";
            //std::cout << "Вторая T0: " << deriv_state[1][5] << "\n";

            // Перевод из кеплеровых координат в декартовы
            init_star_state(x, cur_val[full_size-1], previous_t);

            // Инициализация авто производных
            init_deriv(x);

            // Без массы
            for(int j=0; j< deriv_arr_size-2; j++)
            {
                init_star_state(deriv_state[j], cur_val[full_size-1], previous_t);
                init_deriv(deriv_state[j]);
            }

            //std::cout << "Первая траектория:\n";
            //for(int j = 0; j<6; j++)
            //    printf("%.4e ", deriv_state[0][j]);
            //printf("\n\n");
            //std::cout << "Вторая траектория:\n";
            //for(int j = 0; j<6; j++)
            //    printf("%.4e ", deriv_state[1][j]);
            //printf("\n\n");


            // С массой

            init_star_state(deriv_state[deriv_arr_size-2], cur_val[full_size-1]+(cur_val[full_size-1]/EPS), previous_t);
            init_deriv(deriv_state[deriv_arr_size-2]);

            init_star_state(deriv_state[deriv_arr_size-1], cur_val[full_size-1]-(cur_val[full_size-1]/EPS), previous_t);
            init_deriv(deriv_state[deriv_arr_size-1]);


            while (fgets(buffer, sizeof(buffer), files[file_number]) != NULL)
            {

            // Чтение данных из файла
            sscanf(buffer, "%lf %lf %lf %lf %lf",
                   &t, &ra, &dec, &ra_err, &dec_err);

            // Численное интегирование:
            // Вектор системы
            wrap_integration(x, (t-previous_t)*365.*86400., cur_val[full_size-1], rk_4);

            // Векторы дополнительных орбит для производных
            for (int j=0; j< deriv_arr_size-2; j++)
                wrap_integration(deriv_state[j], (t-previous_t)*365.*86400., cur_val[full_size-1], rk_4);

            // Интегрирование векторов производных по массе (требует eps в wrap_integration)
            wrap_integration(deriv_state[deriv_arr_size-2], (t-previous_t)*365.*86400., cur_val[full_size-1]+(cur_val[full_size-1]/EPS), rk_4);
            wrap_integration(deriv_state[deriv_arr_size-1], (t-previous_t)*365.*86400., cur_val[full_size-1]-(cur_val[full_size-1]/EPS), rk_4);

            // Массив производных
            double deriv[deriv_arr_size];

            // Авто-производные по декартовым элементам (отключены)
            /*
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

            deriv[12] = c/d * x[7];
            deriv[13] = c/d * x[6];
            */

            // Вычисление невязки и проивзодных

            // Невязка
            r_i[0] = c/d* x[1] - ra; // ra под 1
            r_i[1] = c/d* x[0] - dec; // dec под 0

            // Взвешенная сумма квадратов невязок
            sum += pow(r_i[0], 2) / pow(ra_err, 2);
            sum += pow(r_i[1], 2) / pow(dec_err, 2);
            /*
            tmp = 0;

            // Производные (Кроме производной по массе)
            for(int j=0; j<6; j++)
            {
                if (conditions[j] == 1)
                {
                    // Центральные разности
                    // По ra
                    deriv[tmp*2] =
                    c/d*(deriv_state[tmp*2][1] - deriv_state[tmp*2+1][1])/(2*std::abs(cur_val[tmp+size*file_number]/EPS));
                    // по dec
                    deriv[tmp*2+1]=
                    c/d*(deriv_state[tmp*2][0] - deriv_state[tmp*2+1][0])/(2*std::abs(cur_val[tmp+size*file_number]/EPS));
                    // Производные по z и v_z не нужно маштабировать
                    if (j==2 && j==5 && false)
                    {
                        deriv[tmp*2] /= (c/d);
                        deriv[tmp*2+1] /= (c/d);
                    }
                    tmp++;
                }
            }

            //Производные по массе
            deriv[deriv_arr_size-2] =
            c/d*(deriv_state[deriv_arr_size-2][1] - deriv_state[deriv_arr_size-1][1])/(2*std::abs(cur_val[full_size-1]/EPS));
            deriv[deriv_arr_size-1] =
            c/d*(deriv_state[deriv_arr_size-2][0] - deriv_state[deriv_arr_size-1][0])/(2*std::abs(cur_val[full_size-1]/EPS));
            */
            if (true) {
                init_star(deriv_in, file_number);
                int tmp2 = 0;
                for(int j=0; j<6; j++)
                {
                    if (conditions[j] == 1)
                    {
                        deriv_in[j] = cur_val[tmp2 + size*file_number];
                        tmp2++;
                    }
                }
                deriv_in[6] = cur_val[full_size - 1];

                full_analytic(deriv_a, deriv_in, t);
                tmp2 = 0;
                // Аналитические производные
                for(int j=0; j<6; j++)
                {
                    if (conditions[j] == 1)
                    {

                        printf("j = %d. \n RA: %.4e | %.4e \n DEC: %.4e | %.4e\n\n", j, deriv[tmp*2], deriv_a[j*2], deriv[tmp*2+1], deriv_a[j*2+1]);
                        deriv[tmp*2]   = deriv_a[j*2];
                        deriv[tmp*2+1] = deriv_a[j*2+1];
                        tmp2++;
                    }
                }
                deriv[deriv_arr_size-2] = deriv_a[12];
                deriv[deriv_arr_size-1] = deriv_a[13];
            }

            // Техническая переменная для заполнения AtWr
            // Работает как второй счетчик цикла, который срабатывает
            // Только на тех значениях где conditions[j] == 1
            tmp = 0;

            // Заполнение AtWr(betha)
            for(int j=0; j<6; j++)
            {
                if (conditions[j] == 1)
                {
                    AtWr[size*file_number+tmp] +=
                    (1.0/pow(ra_err, 2))*r_i[0]*deriv[tmp*2] + (1.0/pow(dec_err, 2))*r_i[1]*deriv[tmp*2+1];
                    tmp++;
                }

            }

            AtWr[full_size-1] +=
            (1.0/pow(ra_err, 2))*r_i[0]*deriv[deriv_arr_size-2] + (1.0/pow(dec_err, 2))*r_i[1]*deriv[deriv_arr_size-1];


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
                            (1.0/pow(ra_err, 2)) * deriv[t1*2] * deriv[t2*2] +
                            (1.0/pow(dec_err, 2)) * deriv[t1*2+1] * deriv[t2*2+1];
                            t2++;
                        }
                    }

                    // Заполнение правого столбца
                    AtWA[t1 + size*file_number][full_size-1] +=
                    (1.0/pow(ra_err, 2)) * deriv[deriv_arr_size-2] * deriv[t1*2] +
                    (1.0/pow(dec_err, 2)) * deriv[deriv_arr_size-1] * deriv[t1*2+1];

                    // Заполнение нижней строки (Симметрия)
                    AtWA[full_size-1][t1 + size*file_number] = AtWA[t1 + size*file_number][full_size-1];

                    t1++;
                }
            }

            // Заполнение правого нижнего угла
            AtWA[full_size-1][full_size-1] +=
            1.0/pow(ra_err, 2) * deriv[deriv_arr_size-2] * deriv[deriv_arr_size-2] +
            1.0/pow(dec_err, 2) * deriv[deriv_arr_size-1] * deriv[deriv_arr_size-1];


            previous_t = t;

            }
            int t3=0;

            for(int i=0; i<6; i++)
            {
                if (conditions[i] == 1)
                {
                    Preconditioner[size*file_number+t3] = 1./sqrt(AtWA[t3+ size*file_number][t3+ size*file_number]);
                    t3++;
                }
            }
        }



        // Заполнение значений предобуславливателя (масса)

        Preconditioner[full_size-1] = 1./sqrt(AtWA[full_size-1][full_size-1]);

        // Диагональное предобуславливание Якоби

        for(int i=0; i<full_size; i++)
        {
            AtWr[i] *= Preconditioner[i];
            for(int j=0; j<full_size; j++)
                AtWA[i][j] *= Preconditioner[i]*Preconditioner[j];
        }



        std::cout << "------------ITERATION " << i << " -----------------" << std::endl;
        std::cout << std::endl;

        std::cout << "ERROR SUM: " << sum << std::endl;
        std::cout << std::endl;

        if (std::isnan(sum))
            return;

        for(int j=0; j<full_size; j++)
            printf("%.8e ", cur_val[j]);
        std::cout << std::endl;



        std::cout << std::endl;
        std::cout << std::endl;

        for(int j=0; j<full_size; j++)
            printf("%.4e ", AtWr[j]);
        std::cout << std::endl;

        std::cout << std::endl;
        std::cout << std::endl;

        for(int j=0; j<full_size; j++)
        {
            for(int k=0; k< full_size; k++)
                printf("%.4e ", AtWA[j][k]);
            std::cout << std::endl;
        }

        std::cout << std::endl;
        std::cout << std::endl;

        double w[full_size];


        for(int j=0; j<full_size; j++)
        {
            if (j == full_size - 2 || true){
               //AtWA[j][j] *= 2.1;
            }
        }

        solve_eq(AtWA, AtWr, full_size, w);



        for(int j=0; j<full_size; j++) cur_val[j] = cur_val[j] - Preconditioner[j]*w[j];

    }

    // Освобождение памяти
    rk4Free(&rk_4);

    for (int j=0; j<size; j++)
        free(AtWA[j]);
    free(AtWA);

    for (int j=0; j<deriv_arr_size; j++)
        free(deriv_state[j]);

    fclose(files[0]);
    fclose(files[1]);
    fclose(files[2]);
    free(deriv_a);
}
