#!/bin/bash

for i in 500nm 100nm 50nm
do
  echo "Running simulation: $i"
  ./acme2_cyl_par 0 10 1e3 config_shrek/acme2_cyl_par_axonbundle_shrek_square10e6_volumeCorrected_equi_${i}.config &> log_shrek/acme2_cyl_par_axonbundle_shrek_square10e6_equi_${i}.log

  echo "Copying files..."
  cp /import/m1/jpods/dune-ax1/axonbundle_square10e6_equi/config/acme2_cyl_par_axonbundle_volumeCorrected_${i}_equilibrium_basegrid.dgf simulation_states/equilibrium/cyl/yasp_square10mm_axonbundle_volumeCorrected_${i}.dgf

  cp simulation_states/hackepeter/acme2_cyl_par_axonbundle_volumeCorrected_${i}_p0.dat simulation_states/equilibrium/cyl/yasp_square10mm_axonbundle_volumeCorrected_${i}_p0.dat
  echo "...done!"

done


