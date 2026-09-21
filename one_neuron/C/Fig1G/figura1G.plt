reset
set terminal png size 3000,2000 enhanced font 'Times,60'
set output 'Figura1_G.png'
set multiplot layout 2,2
set xrange [0:0.3]
set yrange [0.02:0.06]
set xtics 0.1
set ytics 0.01 offset 0.4
set ylabel "{g_{M}} " offset 1.8 font 'Times,90'
set xlabel "{g_{L}}" offset 0,0.4 font 'Times,90'
########## 1:1
set cbrange [4:16.0]
set cbtics 4 offset -0.8
plot "HH_Freq_RS_I200pA_gT0.400e-3_gL_gM.dat" u 1:2:13 with image t"
########## 2:1
set cbrange [0.3:1.2]
set cbtics 0.3
plot "HH_Freq_RS_I200pA_gT0.400e-3_gL_gM.dat" u 1:2:18 with image t"
########## 2:1
set ylabel "{g_{T}} " offset 1.8 font 'Times,90'
set yrange [0:0.6]
set ytics (0 "0", 0.21 "0.21", 0.4 "0.40", 0.6 "0.6")
set cbrange [4:16.0]
set cbtics 4
set cblabel "{F (Hz)}" offset -6.0,5.25 rotate by 0 left font 'Times,90'
plot "HH_Freq_RS_I200pA_gM0.030e-3_gL_gT.dat" u 1:2:13 with image t"
########## 2:2
set cbrange [0.3:1.2]
set cbtics 0.3
set cblabel "{CV}" offset -6.0,5.25 rotate by 0 left font 'Times,90'
plot "HH_Freq_RS_I200pA_gM0.030e-3_gL_gT.dat" u 1:2:18 with image t"
unset multiplot
