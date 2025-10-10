set terminal wxt size 1280,720 enhanced persist
set xlabel "ΔR.A. (arcsec)"
set ylabel "ΔDec. (arcsec)"

set xrange [*:*] reverse
set xtics 0.1
set ytics 0.05

plot 's2_angles.txt' using 1:2 with points pt 7 ps 0.5 lc rgb "blue" title "S2", \
     's38_angles.txt' using 1:2 with points pt 7 ps 0.5 lc rgb "green" title "S38", \
     's55_angles.txt' using 1:2 with points pt 7 ps 0.5 lc rgb "orange" title "S55/S0-102"

set key top right