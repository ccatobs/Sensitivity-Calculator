#To use this, you MUST be in the Excel Spreadsheet folder and have the am code (I used version 12.0) installed
#This script takes about 15 minutes on my machine; it produces console outputs via am code's automatic outputs as a sign it's running
#! usr/bin/bash

#Values taken from https://arxiv.org/pdf/2007.04262.pdf, feel free to update
CerroChajnantorPWVQ1=0.36
CerroChajnantorPWVQ2=0.67
CerroChajnantorPWVQ3=1.28

#Values approximated by taking Chajnantor Plateau values from https://arxiv.org/pdf/2007.04262.pdf
ActPWVQ1=0.60
ActPWVQ2=1.05
ActPWVQ3=1.98

CCATQ1=$(bc <<<"scale=4; $CerroChajnantorPWVQ1/$ActPWVQ1")
CCATQ2=$(bc <<<"scale=4; $CerroChajnantorPWVQ2/$ActPWVQ2")
CCATQ3=$(bc <<<"scale=4; $CerroChajnantorPWVQ3/$ActPWVQ3")

echo $CCATQ1
echo $CCATQ2
echo $CCATQ3

#Annual
for i in {15..75}
do
    #ACT calculations
    for p in 25 50 75
    do
        am data/ACT_annual_$((p)).amc  0 GHz  1000 GHz  10 MHz  $((i)) deg  1.0 >data/$((p))/ACT_annual_$((p)).$((i)).out
    done

    #Approximations
    am data/ACT_annual_50.amc  0 GHz  1000 GHz  10 MHz  $((i)) deg  0.47058824 >data/25_approx/ACT_annual_25_approximation.$((i)).out
    am data/ACT_annual_50.amc  0 GHz  1000 GHz  10 MHz  $((i)) deg  2.16470588 >data/75_approx/ACT_annual_75_approximation.$((i)).out

    #Steve's approximations
    am data/ACT_annual_50.amc  0 GHz  1000 GHz  10 MHz  $((i)) deg  0.4 >data/Steve/25/ACT_annual_25.$((i)).out
    am data/ACT_annual_50.amc  0 GHz  1000 GHz  10 MHz  $((i)) deg  0.85 >data/Steve/50/ACT_annual_50.$((i)).out
    am data/ACT_annual_50.amc  0 GHz  1000 GHz  10 MHz  $((i)) deg  1.84 >data/Steve/75/ACT_annual_75.$((i)).out

    #CCAT Calculations
    am data/ACT_annual_25.amc  0 GHz  1000 GHz  10 MHz  $((i)) deg  $CCATQ1 >data/CCAT/25/ACT_annual_25.$((i)).out
    am data/ACT_annual_50.amc  0 GHz  1000 GHz  10 MHz  $((i)) deg  $CCATQ2 >data/CCAT/50/ACT_annual_50.$((i)).out
    am data/ACT_annual_75.amc  0 GHz  1000 GHz  10 MHz  $((i)) deg  $CCATQ3 >data/CCAT/75/ACT_annual_75.$((i)).out
done