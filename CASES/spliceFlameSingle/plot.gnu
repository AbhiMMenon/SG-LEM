set yrange [100:2000]
set xrange [0:0.006]
set term pdfcairo #size 640, 480  font "Helvetica,8"
set output 'flame.pdf'


file1 = "results/2.dat"
file2 = "results/100.dat"
file3 = "results/104.dat"
set y2range[-500:5000]
set ylabel "Temperature (K)"
set xlabel "x"
#set y2tics
#set y2label "Mass frac"
plot file1 u 2:8  w lp pi -3 pt 7 ps 0.2   axes x1y1 notitle , \
     file2 u 2:8  w lp pi -3 pt 7 ps 0.2   axes x1y1 notitle , \
     file3 u 2:8  w lp pi -3 pt 7 ps 0.2   axes x1y1 notitle , \
    #file1 u 2:11  w l lt 2 lw 1  axes x1y2 title , \
    #file2 u 2:11  w l lt 2 lw 1  axes x1y2 title , \
    #file3 u 2:11  w l lt 2 lw 1  axes x1y2 title , \
