#!/bin/bash

for i in 1000nm 500nm 100nm 50nm
do
  echo "Running simulation: $i"

  mkdir /import/m1/jpods/dune-ax1/axonbundle_square10e6_volumeCorrected_${i}
  mpirun -bind-to-core -report-bindings -np 10 ./acme2_cyl_par 0 1 20e3 config_shrek/acme2_cyl_par_axonbundle_shrek_square10e6_volumeCorrected_${i}.config &> log_shrek/acme2_cyl_par_axonbundle_shrek_square10e6_volumeCorrected_${i}.log

  echo "...done!"

done


