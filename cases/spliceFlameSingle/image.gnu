#file = 'scalarCInit.dat'
#file2 = 'scalarCInitZ66.dat'
#
#set y2tics
#coPos = 7
#plot file u 1:6 w l axes x1y1 lt 1, file2 u 1:6 w l axes x1y1 lt 3 title "old Z66", \
##file u 1:coPos w l axes x1y2 title "2 step CO", file2 u 1:21 w l axes x1y2 title "Z66 CO"





















set term pdfcairo enhanced dashed size 4,3 font "Times,14"
#set  size 4,3 
set font "Times,11"

#--- --- --- 
file = "scalarReC.dat"
#file2 = "../instanceInit/scalarReC.dat"

pos = 5 
name = 'dcdt'
set output 'map_mono_'.name.'.pdf
load 'default.plt'

set rmargin at screen 0.75
set tmargin at screen 0.8
set lmargin at screen 0.1
set tmargin at screen 0.98
set bmargin at screen 0.2
set ylabel "c"
set xlabel "Re_t"
#set cblabel "dc/dt [s^{-1}]" rotate offset 2,0
set cblabel name rotate offset 2,0
set pm3d noborder map
splot file  u 1:2:pos   notitle w lines palette, file  u 1:2:pos  notitle w pm3d

set term pdfcairo enhanced size 3,3 font "Times,14"
set lmargin at screen 0.2
set output 'graph_mono_'.name.'.pdf
set xlabel "c" offset -2,0
set ylabel name offset 1,0
#set ylabel "dc/dt [s^{-1}]" 
set grid
# set palette model HSV functions (1-gray)*(h2-h1)+h1,1,0.68
set cblabel "Re_t" rotate offset 1,0
#  set cbrange [0:101]
plot for [j=0:11:1] file u 2:pos:1 every :::j::j  w l lw 2 lc palette  notitle, \
#    for [j=0:100:1] file2 u 2:pos:1 every :::j::j  w l lw 2 lc palette  notitle
set grid

