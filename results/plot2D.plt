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
!cat na.dat | grep '#' | wc -l

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

# Write all plots into one file for now, one plot = one page
set output "test.pdf"

set output "flux.pdf"
set title "Membrane flux over time"
set xlabel "Time [ms]"
set ylabel "Membrane flux [mol]/[m^2][s]"
#set output "2d_equi_memb_flux.pdf"
#plot "memb_flux_na.dat" u ($1/1000):2 title "[Na] membrane flux" w lp ls 1,\
#      "memb_flux_k.dat" u ($1/1000):2 title "[K] membrane flux" w lp ls 2
plot "memb_flux_na.dat" u ($1/1000):2 title "[Na] membrane flux" w l ls 1,\
      "memb_flux_k.dat" u ($1/1000):2 title "[K] membrane flux" w l ls 2

set output "flux_sum.pdf"
set title "Sum of membrane fluxes over time"
set xlabel "Time [ms]"
set ylabel "Ratio"      
plot "< join memb_flux_na.dat memb_flux_k.dat" \
  u ($1/1000):($2+$3) title "[Na]+[K] membrane flux" w l ls 1

set output "flux_ratio.pdf"
set title "Ratio of membrane fluxes over time"
plot "< join memb_flux_na.dat memb_flux_k.dat" \
  u ($1/1000):(-$2/$3) title "Ratio [Na]/[K] membrane flux" w l ls 1, \
""  u ($1/1000):(0.13/0.87) title "abs ratio [Na]/[K] leak channels" w l ls 2
           
set output "memb_pot.pdf"
set title "Membrane potential over time"
set xlabel "Time [ms]"
set ylabel "Membrane potential [mV]"                  
#set output "2d_equi_memb_pot.pdf"
#plot "memb_flux_na.dat" u ($1/1000):2 title "[Na] membrane flux" w lp ls 1,\
#      "memb_flux_k.dat" u ($1/1000):2 title "[K] membrane flux" w lp ls 2
plot "memb_pot.dat" u ($1/1000):2 title "Membrane potential" w l ls 3

set output "equi_flux.pdf"
set title "y-component of equilibrium membrane flux over x position"
set xlabel "x position [µm]"
set ylabel "Membrane flux [mol]/[m^2][s]"
plot "< pyextract flux_na.dat x=500" i 100 u 1:3 title "[Na] equilibrium flux y-component" w l ls 1, \
     "< pyextract flux_k.dat x=500" i 100 u 1:3 title "[K] equilibrium flux y-component" w l ls 2, \
     "< pyextract flux_cl.dat x=500" i 100 u 1:3 title "[Cl] equilibrium flux y-component" w l ls 1
     
set output "equi_na.pdf"
set title  "Na concentration equilibrium profile"
set xlabel "x position [µm]"
set ylabel "Na concentration [mM]"
plot "< pyextract na.dat x=500" i 100
      
set output "equi_k.pdf"
set title  "K concentration equilibrium profile"
set xlabel "x position [µm]"
set ylabel "K concentration [mM]"
plot "< pyextract k.dat x=500" i 100

set output "equi_cl.pdf"
set title  "Cl concentration equilibrium profile"
set xlabel "x position [µm]"
set ylabel "Cl concentration [mM]"
plot "< pyextract cl.dat x=500" i 100

set output "equi_cd.pdf"
set title  "Charge density equilibrium profile"
set xlabel "x position [µm]"
set ylabel "Charge density [mM]"
plot "< pyextract cd.dat x=500" i 100

#EOF
