#To use this, you MUST be in the Excel Spreadsheet folder and have the am code (I used version 12.0) installed
#This script takes about 15 minutes on my machine; it produces console outputs via am code's automatic outputs as a sign it's running
#! usr/bin/bash

#Annual
for i in {15..75}
do
    for p in 25 50 75
    do
        am data/ACT_annual_$((p)).amc  0 GHz  1000 GHz  10 MHz  $((i)) deg  1.0 >data/$((p))/ACT_annual_$((p)).$((i)).out
    done
    am data/ACT_annual_50.amc  0 GHz  1000 GHz  10 MHz  $((i)) deg  0.47058824 >data/25_approx/ACT_annual_25_approximation.$((i)).out
    am data/ACT_annual_50.amc  0 GHz  1000 GHz  10 MHz  $((i)) deg  2.16470588 >data/75_approx/ACT_annual_75_approximation.$((i)).out
done