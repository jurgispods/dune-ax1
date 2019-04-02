#!/bin/bash

exit_func() {
  echo "SIGTERM detected, exiting."            
  exit 1
}
trap exit_func SIGTERM SIGINT

dir="/dune/dune-ax1/src"

config="config/square10mm.config"
if [ -z "$1" ]
then
  echo "No config file supplied. Using default file $config"
else
  if [ -f "$dir/$1" ]; then
    config="$dir/$1"
    echo "Using user-supplied config file $config"
  elif [ -f "$dir/config/$1" ]; then
    config="$dir/config/$1"
    echo "Using user-supplied config file $config"
  else
    echo "User-supplied config file $config does not exist! Exiting."
    exit 1
  fi    
fi

cd /dune/dune-ax1/src

echo "Starting simulation using 10 processes."
mpirun --allow-run-as-root -np 10 $dir/acme2_cyl_par 0 1 20e3 $config 
