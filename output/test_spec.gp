reset

#set terminal postscript eps enhanced color dl 3 #mono solid
#font 'Times New Roman, 20'
#set size 1.1, 1.5
#set output "test_prot_spec_1.eps"


set xlabel 'E_{k} [GeV]' font 'Times New Roman, 20' offset 0,-.5
set xrange [5.e-1:2.e3]
#set xrange [1.e0:2.e3]
set logscale x
set format x '10^{%L}'
set xtics in scale 1.5,.75
set xtics font 'Times New Roman, 20' offset 0,-.3
set mxtic 10


set ylabel "E^{2.7}_{k} dN/dE_{k} [GeV^{1.7} m^{-2} s^{-1} sr^{-1}]" font 'Times New Roman, 20' offset -0,0
#set yrange [5e2:2e4] # for proton
set logscale y
set format y '10^{%L}'
set ytics in scale 1.5,.75
set ytics font 'Times New Roman, 20' offset 0,0
set mytic 10



######  "set border" and "unset border" control the display of the graph borders for the plot and splot commands.

set border lw 1.5


######  The margin is the distance between the plot border and the outer edge of the canvas.
######  To alter the distance between the inside of the plot border and the data in the plot itself, see "set offsets".
######  "at screen" indicates that the margin is specified as a fraction of the full drawing area.

#set bmargin at screen
set lmargin at screen .2
#set rmargin at screen
#set tmargin at screen

#set key on left top vertical spacing 2 width 1 height 1 samplen 3  #box #horizontal
set key at screen 0.94,.550 spacing 3 width 1 height 1 samplen 3 # box #horizontal
set key font 'Times New Roman, 10'
set key reverse
set key Left
#unset key


#####    function
mc2    = .938272013
A = 1
Z = 1


Ek(Ri) = sqrt(Ri**2 + mc2**2) -mc2
E(Ri)  = sqrt(Ri**2 + mc2**2)
Ri(Ek) = sqrt((Ek+mc2)**2 -mc2**2)

M2p7(Ek) = (Ek)**(2.7)
E2(x1,x2) = exp((log(x1)+log(x2))/2.)

sigma(x1,x2) = sqrt(x1*x1+x2*x2)

p(Ek) = sqrt((Ek+m)**2-m**2)


#plot '/home/liuwei/Research/PROGRAM/CRPropagation/data/AMS-02/AMS-02_prot_2015.dat' u (Ek(sqrt($1*$2))):(Ek(sqrt($1*$2))**2.7*$3*$10):(Ek(sqrt($1*$2))**2.7*sigma($4,$9)*$10)  t 'AMS-02'  w yerrorbars lw 3 pt 5 lc rgb 'red', \

plot 'example/example.dat' u 1:($1**2.7*$22) w l, \
     'example1/example1.dat' u 1:($1**2.7*$22) w l, \