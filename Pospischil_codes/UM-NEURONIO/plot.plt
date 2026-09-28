reset
set terminal postscript eps size 5,8 enhanced color font 'Helvetica,22'
set output 'Figura_g1_gM0.03_delay0.5.eps'
set multiplot layout 2,1
set xlabel "{g_{ex}}" offset 0,0.3 font 'Helvetica,40';
set ylabel "{I_{ext}} " offset 1.3 font 'Helvetica,40'
set xrange [0.00:1.5]
set yrange [0.135:0.28]
set xtics  0.3
set ytics  0.07 offset 0.5
set cbtics 8.0 offset -1.0
set cbrange [0:40.0]
set cblabel "{FR} " offset 0.0,0.0 rotate by 0 left font 'Helvetica,40'
plot "HH_RS_g1_I_gexc_gM0.03_gL0.0_gT0.0_delay0.5.dat" u 4:2:7 w image t"
set cbtics 0.1 offset -1.0
set cbrange [0:1.5]
set cblabel "{CV} " offset 0.0,0.0 rotate by 0 left font 'Helvetica,40'
plot "HH_RS_g1_I_gexc_gM0.03_gL0.0_gT0.0_delay0.5.dat" u 4:2:6 w image t"
unset multiplot
