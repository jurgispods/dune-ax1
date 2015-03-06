#!/bin/bash

set -e

if [ -e "cond_numbers.dat" ]
then
	rm cond_numbers.dat
fi

for file in $(ls octave/nernst_planck_mat*)
do
	./cond_number.m $file | tee -a cond_numbers.dat
done

