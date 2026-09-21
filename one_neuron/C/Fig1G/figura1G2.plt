reset
set terminal png size 3000,2000 enhanced font 'Times,45'
set output 'Figura1_G2.png'
set multiplot layout 2,2
set xrange [0:0.4]
set yrange [0:0.07]
set xlabel "{g_{L}}" offset 0,0.8 font 'Times,52'
set ylabel "{g_{M}} " offset 0.8 font 'Times,52'
########## 1:1
set cbrange [0:30.0]
set cblabel "{FR (Hz)}" offset 0,0.8 font 'Times,52'
plot "HH_Freq_RS_I300pA_gT0.400e-3_gL_gM.dat" u 1:2:13 with image t"
########## 2:1
set cbrange [0.2:1.5]
set cblabel "{CV}" offset 0,0.8 font 'Times,52'
plot "HH_Freq_RS_I300pA_gT0.400e-3_gL_gM.dat" u 1:2:18 with image t"
########## 2:1
set cbrange [0:0.2]
set cblabel "Adpt" offset 0,0.8 font 'Times,52'
plot "HH_Freq_RS_I300pA_gT0.400e-3_gL_gM.dat" u 1:2:21 with image t" 
########## 2:2
set cbrange [24:26.0]
set cblabel "T_1" offset 0,0.8 font 'Times,52'
plot "HH_Freq_RS_I300pA_gT0.400e-3_gL_gM.dat" u 1:2:16 with image t"
unset multiplot
