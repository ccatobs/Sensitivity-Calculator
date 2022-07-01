#To use this, you must be in the Excel Spreadsheet folder and have the am code installed. It takes approximately 2 minutes to run, and produces no output to stdout. The only things created should be files inside folders in data/
#! usr/bin/bash

#Values taken from https://arxiv.org/pdf/2007.04262.pdf
CerroChajnantorPWVQ1=0.36
CerroChajnantorPWVQ2=0.67
CerroChajnantorPWVQ3=1.28

#PWV values in ACT configuration files
ACTConfigQ1=$(python3 data/pwvCalculator.py 25)
ACTConfigQ2=$(python3 data/pwvCalculator.py 50)
ACTConfigQ3=$(python3 data/pwvCalculator.py 75)

#Ratio of PWVs used in CCAT site data generation
CerroConfigQ1=$(bc <<<"scale=10; $CerroChajnantorPWVQ1/$ACTConfigQ1")
CerroConfigQ2=$(bc <<<"scale=10; $CerroChajnantorPWVQ2/$ACTConfigQ2")
CerroConfigQ3=$(bc <<<"scale=10; $CerroChajnantorPWVQ3/$ACTConfigQ3")

#Create base files so ACT's PWVs are known
for p in 25 50 75
do
    am data/ACT_annual_$((p)).amc  0 GHz  1000 GHz  10 MHz  45 deg  1.0 >data/$((p))/ACT_annual_$((p)).45.out 2>data/$((p))/ACT_annual_$((p)).45.err
done

for i in {15..75}
do
    #CerroConfig Calculations
    am data/ACT_annual_25.amc  0 GHz  1000 GHz  10 MHz  $((i)) deg  $CerroConfigQ1 >data/CerroConfig/25/ACT_annual_25.$((i)).out 2>data/CerroConfig/25/ACT_annual_25.$((i)).err
    am data/ACT_annual_50.amc  0 GHz  1000 GHz  10 MHz  $((i)) deg  $CerroConfigQ2 >data/CerroConfig/50/ACT_annual_50.$((i)).out 2>data/CerroConfig/50/ACT_annual_50.$((i)).err
    am data/ACT_annual_75.amc  0 GHz  1000 GHz  10 MHz  $((i)) deg  $CerroConfigQ3 >data/CerroConfig/75/ACT_annual_75.$((i)).out 2>data/CerroConfig/75/ACT_annual_75.$((i)).err
done