#!/bin/bash

vars="minimalMembraneDistance"

for file in $( ls *.config );
do
  echo "Processing file $file..."
  
  for var in $vars;
  do
    echo " - refactoring variable $var..."
    
    # sed magic happening here
    sed -i "s/\(^${var}\)\s*=\s*\([^']*\)\s*$/\1 = 709/g" $file
    sed -i "/\(#${var}\)\s*=\s*\(1214\)\s*$/d" $file
  done
  
done

