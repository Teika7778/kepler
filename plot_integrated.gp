set terminal wxt size 720,720 enhanced persist
set xlabel "x"
set ylabel "y"

set xrange [*:*] reverse
set xtics 0.1
set ytics 0.05

plot 'integration_s2.txt' using 1:2 with points pt 7 ps 0.5 lc rgb "blue" title "S2", \
     'integration_s38.txt' using 1:2 with points pt 7 ps 0.5 lc rgb "green" title "S38", \
     'integration_s55.txt' using 1:2 with points pt 7 ps 0.5 lc rgb "orange" title "S55/S0-102"

set key top right