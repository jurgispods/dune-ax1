#!/usr/bin/gnuplot -persist
set terminal pdfcairo \
  enhanced color \
  dashed \
  font "Gill Sans,8" \
  linewidth 4 \
  rounded #size 10cm,8cm
  
#set style line 
# set output 'output.ps'
#set xlabel "x" 0.000000,0.000000  ""
#set ylabel "y=exp(-x)" 0.000000,0.000000  ""
#set title "Pade approximation" 0.000000,0.000000  ""
#set xrange [ 0 : 2 ] noreverse nowriteback
#set yrange [ 0 : 1 ] noreverse nowriteback
#set mxtics 5.000000
#set mytics 5.000000
#set xtics border mirror norotate 1
#set ytics border mirror norotate 0.5
#!cat na.dat | grep '#' | wc -l

# Line style for axes
set style line 80 lt rgb "#808080"

# Line style for grid
set style line 81 lt 0  # dashed
set style line 81 lt rgb "#808080"  # grey

set grid back linestyle 81
set border 3 back linestyle 80 # Remove border on top and right.  These
             # borders are useless and make it harder
             # to see plotted lines near the border.
    # Also, put it in grey; no need for so much emphasis on a border.
set xtics nomirror
set ytics nomirror

#set log x
#set mxtics 10    # Makes logscale look good.


# Line styles: try to pick pleasing colors, rather
# than strictly primary colors or hard-to-see colors
# like gnuplot's default yellow.  Make the lines thick
# so they're easy to see in small plots in papers.
set style line 1 lt rgb "#A00000" lw 2 pt 1
set style line 2 lt rgb "#00A000" lw 2 pt 6
set style line 3 lt rgb "#5060D0" lw 2 pt 2
set style line 4 lt rgb "#F25900" lw 2 pt 9

set style data lines
set key top right

set output "isolines.pdf"

reset
#set xrange [-5:5]
#set yrange [-5:5]
set isosample 250, 250
set table 'test.dat'
splot "pot.dat" i 200
unset table


set contour base
set cntrparam level incremental -70, 5, 60
unset surface
set table 'cont.dat'
splot "pot.dat" i 200
unset table

reset
#set xrange [-5:5]
#set yrange [-5:5]
unset key
set palette rgbformulae 33,13,10
p 'test.dat' with image, 'cont.dat' w l lt -1 lw 1.5

