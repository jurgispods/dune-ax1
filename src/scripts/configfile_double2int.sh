#!/bin/bash

vars="checkpointInterval n_dy_min n_dy_cell_min"

for file in $( ls *.config );
do
  echo "Processing file $file..."
  
  for var in $vars;
  do
    echo " - refactoring variable $var..."
    sed -i "s/\(${var}\)\s*=\s*\([^']*\)\.\s*$/\1 = \2/g" $file
  done
  
done

