reset
set terminal png size 5500,2000 enhanced font 'Times,60'
set output 'FiguraS2.png'
set multiplot layout 2,4
set xrange [0:0.4]
set yrange [0.01:0.06]
set xtics 0.1
set ytics 0.02 offset 0.4
set ylabel "{g_{M}} " offset 1.8 font 'Times,90'
set xlabel "{g_{L}}" offset 0,0.4 font 'Times,90'
########## 1:1
set ylabel "{g_{T}} " offset 1.8 font 'Times,90'
set yrange [0:0.6]
set ytics (0 "0", 0.21 "0.21", 0.4 "0.40", 0.6 "0.6")
set cbrange [10:70.0]
set cbtics 20
set cblabel "{F (Hz)}" offset -6.0,5.25 rotate by 0 left font 'Times,90'
plot "HH_Freq_RS_I200pA_gM0.030e-3_gL_gT.dat" u 1:2:13 with image t"
########## 2:2
set cbrange [0.0:3.0]
set cbtics 1.0
set cblabel "{CV}" offset -6.0,5.25 rotate by 0 left font 'Times,90'
plot "HH_Freq_RS_I200pA_gM0.030e-3_gL_gT.dat" u 1:2:18 with image t"
########## 2:3
set cbrange [0.0:0.2]
set cbtics 0.1
set cblabel "{A_{10}}" offset -6.0,5.25 rotate by 0 left font 'Times,90'
plot "HH_Freq_RS_I200pA_gM0.030e-3_gL_gT.dat" u 1:2:21 with image t"
########## 2:4
set cbrange [0.0:0.1]
set cbtics 0.05
set cblabel "{A_{20}}" offset -6.0,5.25 rotate by 0 left font 'Times,90'
plot "HH_Freq_RS_I200pA_gM0.030e-3_gL_gT.dat" u 1:2:23 with image t"
########## 2:1
set ylabel "{g_{T}} " offset 1.8 font 'Times,90'
set yrange [0:0.6]
set ytics (0 "0", 0.21 "0.21", 0.4 "0.40", 0.6 "0.6")
set cbrange [10:70.0]
set cbtics 20
set cblabel "{F (Hz)}" offset -6.0,5.25 rotate by 0 left font 'Times,90'
plot "HH_Freq_RS_I200pA_gM0.050e-3_gL_gT.dat" u 1:2:13 with image t"
########## 2:2
set cbrange [0.0:3.0]
set cbtics 1.0
set cblabel "{CV}" offset -6.0,5.25 rotate by 0 left font 'Times,90'
plot "HH_Freq_RS_I200pA_gM0.050e-3_gL_gT.dat" u 1:2:18 with image t"
########## 2:3
set cbrange [0.0:0.2]
set cbtics 0.1
set cblabel "{A_{10}}" offset -6.0,5.25 rotate by 0 left font 'Times,90'
plot "HH_Freq_RS_I200pA_gM0.050e-3_gL_gT.dat" u 1:2:21 with image t"
########## 2:4
set cbrange [0.0:0.1]
set cbtics 0.05
set cblabel "{A_{20}}" offset -6.0,5.25 rotate by 0 left font 'Times,90'
plot "HH_Freq_RS_I200pA_gM0.050e-3_gL_gT.dat" u 1:2:23 with image t"
unset multiplot
