#!/bin/bash

exit_func() {
  echo "SIGTERM detected, exiting."            
  exit 1
}
trap exit_func SIGTERM SIGINT

dir="/dune/dune-ax1/src"

echo "args: $@"

if [ "$#" -ne 0 ] && [ "$#" -gt 6 ]; then
    echo "Usage: docker run dune-ax1 <executable> <refinement_level> <dt_start> <t_end> <config> <n_proc>"
    echo " - executable: 'acme2_cyl_par' (unmyelinated setup) or 'acme2_cyl_par_myelin' (myelinated setup), default: acme2_cyl_par"
    echo " - refinement_level: Number of global grid refinements, default: 0"
    echo " - dt_start: Initial time step (in microseconds), default: 0"
    echo " - t_end: Simulation end time (in microseconds), default: 20e3"
    echo " - config: Config file to use, default: config/square10mm.config"
    echo " - n_proc: Number of processors to use for parallel run, default: 10"
fi

executable="acme2_cyl_par"
if [ -z "$1" ]
then
  echo "Executable not specified. Using default $executable"
else
  executable="$1"
fi

refinement_level="0"
if [ -z "$2" ]
then
  echo "Refinement level not specified. Using default $refinement_level"
else
  refinement_level="$2"
fi

dt_start="1"
if [ -z "$3" ]
then
  echo "Initial time step not specified. Using default $dt_start"
else
  dt_start="$3"
fi

t_end="20e3"
if [ -z "$4" ]
then
  echo "End time not specified. Using default $t_end"
else
  t_end="$4"
fi

config="config/square10mm.config"
if [ -z "$5" ]
then
  echo "No config file supplied. Using default file $config"
else
  if [ -f "$dir/$5" ]; then
    config="$dir/$5"
    echo "Using user-supplied config file $config"
  elif [ -f "$dir/config/$5" ]; then
    config="$dir/config/$5"
    echo "Using user-supplied config file $config"
  else
    echo "User-supplied config file $config does not exist! Exiting."
    exit 1
  fi    
fi

n_proc="10"
if [ -z "$6" ]
then
  echo "Number of processors not specified. Using default $n_proc"
else
  n_proc="$6"
fi


cd /dune/dune-ax1/src

echo "Starting simulation using $n_proc processes."
mpirun --allow-run-as-root -np $n_proc $executable $refinement_level $dt_start $t_end $config
