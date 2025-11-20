# Непосредственное вычисление разницы
set terminal pngcairo size 1200,800 enhanced font 'Arial,12'
set output "diriv_difference.png"

set multiplot layout 2,1 title "Модуль разницы численной и аналитической производными в рамках одной орбиты"

set title "Разница ra"
set xlabel "Номер точки"
set grid

plot "< paste num_deriv_value.txt analytic_deriv_value.txt" using 0:(abs($1-$3)) with linespoints title 'Разница ra'

set title "Разница dec"
set xlabel "Номер точки"

plot "< paste num_deriv_value.txt analytic_deriv_value.txt" using 0:(abs($2-$4)) with linespoints title 'Разница dec'

unset multiplot