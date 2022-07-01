#To use this, you must be in the Excel Spreadsheet folder and have the am code (I used version 12.0) installed
#This script takes about 2 minutes on my machine, so it's probably not broken
#! usr/bin/bash

#Values taken from https://arxiv.org/pdf/2007.04262.pdf
CerroChajnantorPWVQ1=0.36
CerroChajnantorPWVQ2=0.67
CerroChajnantorPWVQ3=1.28

#PWV values in ACT configuration files, collected by data/pwvCalculator.py
ACTConfigQ1=0.457488
ACTConfigQ2=0.9315249999999999
ACTConfigQ3=2.23472

CerroConfigQ1=$(bc <<<"scale=10; $CerroChajnantorPWVQ1/$ACTConfigQ1")
CerroConfigQ2=$(bc <<<"scale=10; $CerroChajnantorPWVQ2/$ACTConfigQ2")
CerroConfigQ3=$(bc <<<"scale=10; $CerroChajnantorPWVQ3/$ACTConfigQ3")

#Get PWV values from config files
for p in 25 50 75
do
    am data/ACT_annual_$((p)).amc  0 GHz  1000 GHz  10 MHz  45 deg  1.0 >data/$((p))/ACT_annual_$((p)).45.out 2>data/$((p))/ACT_annual_$((p)).45.err
done

#Annual
for i in {15..75}
do
    #CerroConfig Calculations
    am data/ACT_annual_25.amc  0 GHz  1000 GHz  10 MHz  $((i)) deg  $CerroConfigQ1 >data/CerroConfig/25/ACT_annual_25.$((i)).out 2>data/CerroConfig/25/ACT_annual_25.$((i)).err
    am data/ACT_annual_50.amc  0 GHz  1000 GHz  10 MHz  $((i)) deg  $CerroConfigQ2 >data/CerroConfig/50/ACT_annual_50.$((i)).out 2>data/CerroConfig/50/ACT_annual_50.$((i)).err
    am data/ACT_annual_75.amc  0 GHz  1000 GHz  10 MHz  $((i)) deg  $CerroConfigQ3 >data/CerroConfig/75/ACT_annual_75.$((i)).out 2>data/CerroConfig/75/ACT_annual_75.$((i)).err
done