#To use this, you MUST be in the data folder and have the am code (I used version 12.0) installed
#! usr/bin/bash

for i in {0..89}
do
    for p in 25 50 75
    do
        var=`am ACT_annual_$((p)).amc  0 GHz  1000 GHz  10 MHz  $((i)) deg  1.0 >ACT_annual_$((p)).$((i)).out`
        echo $var
    done

    var=`am ACT_annual_50.amc  0 GHz  1000 GHz  10 MHz  $((i)) deg  0.47058824 >ACT_annual_25_approximation.$((i)).out`
    echo $var

    var=`am ACT_annual_50.amc  0 GHz  1000 GHz  10 MHz  $((i)) deg  2.16470588 >ACT_annual_75_approximation.$((i)).out`
    echo $var
done
