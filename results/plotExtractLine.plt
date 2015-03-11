#!/bin/bash

if test $# -lt 2
then
  echo "Usage: ./plot.plt <gnuplot data file> <node index>"
  exit
fi

arr=$(echo $1 | tr "/" "\n")

folder=""
i=0
for x in $arr
do
  if [[ $i == "0" ]]
  then
    folder=$x
  fi
  echo "> [$x]"
  i=`expr $i + 1`
done

nLines=`wc -l $1`
nTimeSteps=`grep "#" $1 | wc -l`

echo "Number of lines: $nLines"
echo "Number of time steps: $nTimeSteps"

python extractLine.py $1 $2 > temp1.dat
#python extractIndex.py $1 pos=50.0 > temp2.dat
#plotData=`python extractIndex.py $1 $2`
#echo $plotData

gnuplot -persist << EOF
set terminal wxt
#plot "-"
#!python extractLine.py $1 $2
#e
plot "temp1.dat" title "position y=$2" w l
#replot "temp2.dat" using 1:3 title "position @ index 50." w l
#replot "${folder}/memb_pot.dat" using 1:2 title "membrane" w l
EOF

