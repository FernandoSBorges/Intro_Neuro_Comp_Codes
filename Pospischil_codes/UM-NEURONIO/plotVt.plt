reset
set terminal png size 800,800 font 'Times,22'
set output 'Figura_g1_gM0.03_delay0.5.png'
set multiplot layout 2,1
set ylabel "V (mV) " offset 1.3  font 'Times,25'
set xrange [0.0:1000.0]
set yrange [-80.0:40.0]
set xtics  500.0
set ytics  20.0 offset 0.5
plot "HH_N1.dat" u 1:2 with lines lw 2.5 t"
set ylabel "m" offset 1.3   font 'Times,28'
set xlabel "time (ms)" offset 0,0.3 font 'Times,25'
set yrange [0.0:0.4]
set ytics  0.2 offset 0.5
plot "HH_N1.dat" u 1:3 with lines lw 2.5 t"
unset multiplot
