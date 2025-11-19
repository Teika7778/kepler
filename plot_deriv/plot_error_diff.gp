set terminal pngcairo size 1200,600 enhanced font 'Arial,12'
set output "error_sum_difference.png"

set title "Сравнение суммы ошибок по итерациям"
set xlabel "Номер итерации"
set ylabel "Сумма ошибок"
set grid

plot "num_deriv_error_sum.txt" with linespoints title 'Численные производные', \
     "analytic_deriv_error_sum.txt" with linespoints title 'Аналитические производные'

plot "< paste num_deriv_error_sum.txt analytic_deriv_error_sum.txt" using 1:($2-$4) with linespoints title 'Разница'