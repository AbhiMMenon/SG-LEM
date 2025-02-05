#set term gif animate size 640, 480  font "Helvetica,8"
#set output 'testAnimate.gif'


tStep = 1E-4
nStep = 1200
#--- --- --- 

#set grid
#set xrange[0:.01]
#set yrange[200:2000]
#set y2range[0:0.03]
set ylabel "Temperature (K)"
set xlabel "x"
set y2tics
set y2label "Mass frac"
do for [i=1:nStep] {
    file = sprintf("./results/%i.dat", i)
    comm = sprintf("head -n 1 %s", file)
    s    = system(comm)
    f = sprintf("Time = %s s", word(s,3))
    set title f
    plot file  u 2:24  w lp lt 1 lw 1  axes x1y1 title "T"
}

#plot for [i=2:200] sprintf("./results/%i.dat",i) u 2:8 w l  notitle

