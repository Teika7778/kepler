#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>

#include "constans.hpp"
#include "struct.hpp"
#include "normalize.hpp"
#include "transform.hpp"
#include "diff.hpp"
#include "gauss_newton.hpp"
#include "integration.hpp"

int main()
{
    // Открываем файл с входными данными
    std::ifstream input_file("user_files/input.txt");
    if (!input_file.is_open()) {
        std::cerr << "Error: cant open file input.txt!" << std::endl;
        return 1;
    }

    int GN_num_iter;
    double alpha;

    // 1. Читаем настройки оптимизатора
    input_file >> GN_num_iter >> alpha;

    // 2. Читаем флаги определяемых параметров (a, e, w, omega, i, T0)
    int conditions[6];
    for (int j = 0; j < 6; j++) {
        input_file >> conditions[j];
    }

    int num_star = 3;
    int full_size = 6 * num_star + 1; // 19 параметров (18 для орбит + 1 масса)
    std::vector<double> parameters(full_size);

    // 3. Читаем параметры звезд (S2, S38, S55)
    for (int j = 0; j < 6 * num_star; j++)
        input_file >> parameters[j];

    // 4. Читаем массу черной дыры
    input_file >> parameters[full_size - 1];

    // 5. Читаем процент ошибки и возмущаем на нее
    for(int j=0; j<7; j++)
    {
        double error;

        input_file >> error;

        if (j == 6)
            parameters[18] *= (1 + error / 100);
        else
        {
            parameters[j] *= (1 + error / 100);
            parameters[j + 6] *= (1 + error / 100);
            parameters[j + 12] *= (1 + error / 100);
        }
        
    }
        

    std::cout << "Data loaded" << std::endl;
    
    input_file.close();

    gauss_newton(parameters.data(), num_star, conditions, GN_num_iter, alpha);

    // --- ВЫВОД РЕЗУЛЬТАТА В result.txt ---
    std::ofstream res_file("user_files/result.txt");
    if (res_file.is_open()) {
        
        // Массив названий звезд для удобства
        std::string star_names[] = {"S2", "S38", "S55"};
        
        for (int s = 0; s < num_star; s++) {
            res_file << "Параметры " << star_names[s] << ":" << std::endl;
            res_file << "a:     " << parameters[s * 6 + 0] << std::endl;
            res_file << "e:     " << parameters[s * 6 + 1] << std::endl;
            res_file << "w:     " << parameters[s * 6 + 2] << std::endl;
            res_file << "omega: " << parameters[s * 6 + 3] << std::endl;
            res_file << "i:     " << parameters[s * 6 + 4] << std::endl;
            res_file << "T0:    " << parameters[s * 6 + 5] << std::endl;
            res_file << "----------------------------" << std::endl;
        }

        res_file << "\nМасса черной дыры:" << std::endl;
        res_file << parameters[full_size - 1] << std::endl;
        
        res_file.close();

        if (std::isnan(parameters[0]))
            std::cout << "Result is nan" << std::endl;
        else
            std::cout << "Result is written to result.txt" << std::endl;
    }

    return 0;
}