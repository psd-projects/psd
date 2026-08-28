set multiplot

set tics font "times,15"
set key font "times,15" bottom
set xlabel 'Time Steps' font "times,15"
set ylabel 'Pressure_{normalized}' font "Helvetica,15"


plot 'result-drucker-prager-mfront.data' w  lp pt 8  ps 3  t 'PSD-mfront', 'result-drucker-prager.data' w  lp pt 6  ps 3  t 'PSD', './ref/ref_druker_prager.data' u 1:2  w lp pt 7 ps 1 t 'Ref'



set tics font "times,10"
set key font "times,10" bottom
set xlabel '' font "times,10"
set ylabel 'Error %' font "Helvetica,10"
set ytic 0.0001 
set size 0.65, 0.4
set origin 0.3, 0.3

plot '< paste result-drucker-prager-mfront.data ./ref/ref_druker_prager.data' u 1:(abs($2-$4)*100/$4) w l lw 3 t ''

unset multiplot

set term postscript
set output "comp-drucker-prager.eps"
load "gnuplot.input"
