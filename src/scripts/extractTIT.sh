#!/bin/bash

file=$1

out1=`echo $file | sed -e "s/log/tit1/"`
out2=`echo $file | sed -e "s/log/tit2/"`

time1=`echo $file | sed -e "s/log/time1/"`
time2=`echo $file | sed -e "s/log/time2/"`
time3=`echo $file | sed -e "s/log/time3/"`
time4=`echo $file | sed -e "s/log/time4/"`
time5=`echo $file | sed -e "s/log/time5/"`
time6=`echo $file | sed -e "s/log/time6/"`
time7=`echo $file | sed -e "s/log/time7/"`
time8=`echo $file | sed -e "s/log/time8/"`
time9=`echo $file | sed -e "s/log/time9/"`

grep -Hn --color "TIT1" $file | sed "s/.*TIT1: //" &> $out1
grep -Hn --color "TIT2" $file | sed "s/.*TIT2: //" &> $out2

grep -Hn --color "time1" $file | sed "s/.*time1: //" &> $time1
grep -Hn --color "time2" $file | sed "s/.*time2: //" &> $time2
grep -Hn --color "time3" $file | sed "s/.*time3: //" &> $time3
grep -Hn --color "time4" $file | sed "s/.*time4: //" &> $time4
grep -Hn --color "time5" $file | sed "s/.*time5: //" &> $time5
grep -Hn --color "time6" $file | sed "s/.*time6: //" &> $time6
grep -Hn --color "time7" $file | sed "s/.*time7: //" &> $time7
grep -Hn --color "time8" $file | sed "s/.*time8: //" &> $time8
grep -Hn --color "time9" $file | sed "s/.*time9: //" &> $time9

