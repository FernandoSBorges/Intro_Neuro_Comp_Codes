reset
set terminal postscript eps size 24,16 enhanced color font 'Times,45'
set output 'Figura-Isyn200_gM0.05_.eps'
set multiplot layout 3,3
set xrange [0:0.4]
set yrange [0:0.8]
########## 1:1
set cbrange [0:0.08]
plot "HH_Freq_RS_I200pA_gM0.050e-3_gL_gT.dat" u 1:2:22 with image t" 
########## 2:1
set cbrange [0:0.08]
plot "HH_Freq_RS_I200pA_gM0.050e-3_gL_gT.dat" u 1:2:23 with image t" 
########## 3:1
set cbrange [0:0.2]
plot "HH_Freq_RS_I200pA_gM0.050e-3_gL_gT.dat" u 1:2:21 with image t" 
########## 1:2
set cbrange [5:15.0]
plot "HH_Freq_RS_I200pA_gM0.050e-3_gL_gT.dat" u 1:2:13 with image t"
########## 2:2
set cbrange [0.2:1.2]
plot "HH_Freq_RS_I200pA_gM0.050e-3_gL_gT.dat" u 1:2:18 with image t"
########## 3:2
set cbrange [0:0.4]
plot "HH_Freq_RS_I200pA_gM0.050e-3_gL_gT.dat" u 1:2:20 with image t" 
########## 1:3
set cbrange [0:0.08]
plot "HH_Freq_RS_I200pA_gM0.050e-3_gL_gT.dat" u 1:2:24 with image t" 
########## 2:3
set cbrange [0:0.08]
plot "HH_Freq_RS_I200pA_gM0.050e-3_gL_gT.dat" u 1:2:25 with image t" 
########## 3:3
set cbrange [0:0.35]
plot "HH_Freq_RS_I200pA_gM0.050e-3_gL_gT.dat" u 1:2:19 with image t" 
unset multiplot
