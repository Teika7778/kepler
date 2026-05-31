set terminal pngcairo size 1200,800 enhanced font 'Arial,12'
set output "derivatives.png"

set multiplot layout 2,1 title "Сравнение производных в рамках одной орбиты"

set title "Значение по ra"
set xlabel "Номер точки"

set grid

plot "num_deriv_value.txt" using 0:1 with lines title 'Численная производная', \
     "analytic_deriv_value.txt" using 0:1 with lines title 'Аналитическая производная'

set title "Значение по dec"
set xlabel "Номер точки"

plot "num_deriv_value.txt" using 0:2 with lines title 'Численная производная', \
     "analytic_deriv_value.txt" using 0:2 with lines title 'Аналитическая производная'

unset multiplot