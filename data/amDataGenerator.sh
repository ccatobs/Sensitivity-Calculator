#! usr/bin/bash

for i in {0..89}
do
    for p in 25 50 75
    do
        var=`am ACT_annual_$((p)).amc  0 GHz  1000 GHz  10 MHz  $((i)) deg  1.0 >ACT_annual_$((p)).$((i)).out`
        echo $var
    done
done

