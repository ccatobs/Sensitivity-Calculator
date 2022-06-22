#To use this, you MUST be in the Excel Spreadsheet folder and have the am code (I used version 12.0) installed
#This script takes about 15 minutes on my machine; it produces console outputs via am code's automatic outputs as a sign it's running
#! usr/bin/bash

#Values taken from https://arxiv.org/pdf/2007.04262.pdf
CerroChajnantorPWVQ1=0.36
CerroChajnantorPWVQ2=0.67
CerroChajnantorPWVQ3=1.28

#ACT PWV values approximated by taking Chajnantor Plateau values from https://arxiv.org/pdf/2007.04262.pdf
ChajnantorPlateauQ1=0.60
ChajnantorPlateauQ2=1.05
ChajnantorPlateauQ3=1.98

#PWV Values from the APEX telescope site, data from https://www.apex-telescope.org/apex-dashboard/d/MQgvc4Onz and calculated by data/pwvCalculator.py
APEXQ1=0.7284166666666666
APEXQ2=1.2816666666666667
APEXQ3=2.28

#PWV values in ACT configuration files, collected by data/pwvCalculator.py
ACTConfigQ1=0.70588868
ACTConfigQ2=1.4705762629999999
ACTConfigQ3=3.55953805

#Calculate PWV ratios between ACT and CCAT sites
CerroPlateauQ1=$(bc <<<"scale=10; $CerroChajnantorPWVQ1/$ChajnantorPlateauQ1")
CerroPlateauQ2=$(bc <<<"scale=10; $CerroChajnantorPWVQ2/$ChajnantorPlateauQ2")
CerroPlateauQ3=$(bc <<<"scale=10; $CerroChajnantorPWVQ3/$ChajnantorPlateauQ3")

CerroAPEXQ1=$(bc <<<"scale=10; $CerroChajnantorPWVQ1/$APEXQ1")
CerroAPEXQ2=$(bc <<<"scale=10; $CerroChajnantorPWVQ2/$APEXQ2")
CerroAPEXQ3=$(bc <<<"scale=10; $CerroChajnantorPWVQ3/$APEXQ3")

CerroConfigQ1=$(bc <<<"scale=10; $CerroChajnantorPWVQ1/$ACTConfigQ1")
CerroConfigQ2=$(bc <<<"scale=10; $CerroChajnantorPWVQ2/$ACTConfigQ2")
CerroConfigQ3=$(bc <<<"scale=10; $CerroChajnantorPWVQ3/$ACTConfigQ3")


#Annual
for i in {15..75}
do
    #ACT calculations
    for p in 25 50 75
    do
        am data/ACT_annual_$((p)).amc  0 GHz  1000 GHz  10 MHz  $((i)) deg  1.0 >data/$((p))/ACT_annual_$((p)).$((i)).out
    done

    #Approximations
    am data/ACT_annual_50.amc  0 GHz  1000 GHz  10 MHz  $((i)) deg  0.47058824 >data/25_approx/ACT_annual_25.$((i)).out
    am data/ACT_annual_50.amc  0 GHz  1000 GHz  10 MHz  $((i)) deg  2.16470588 >data/75_approx/ACT_annual_75.$((i)).out

    #Steve's approximations
    am data/ACT_annual_50.amc  0 GHz  1000 GHz  10 MHz  $((i)) deg  0.51 >data/Steve/25/ACT_annual_25.$((i)).out
    am data/ACT_annual_50.amc  0 GHz  1000 GHz  10 MHz  $((i)) deg  0.95 >data/Steve/50/ACT_annual_50.$((i)).out
    am data/ACT_annual_50.amc  0 GHz  1000 GHz  10 MHz  $((i)) deg  1.81 >data/Steve/75/ACT_annual_75.$((i)).out

    #CerroPlateau Calculations
    am data/ACT_annual_25.amc  0 GHz  1000 GHz  10 MHz  $((i)) deg  $CerroPlateauQ1 >data/CerroPlateau/25/ACT_annual_25.$((i)).out
    am data/ACT_annual_50.amc  0 GHz  1000 GHz  10 MHz  $((i)) deg  $CerroPlateauQ2 >data/CerroPlateau/50/ACT_annual_50.$((i)).out
    am data/ACT_annual_75.amc  0 GHz  1000 GHz  10 MHz  $((i)) deg  $CerroPlateauQ3 >data/CerroPlateau/75/ACT_annual_75.$((i)).out

    #CerroAPEX Calculations
    am data/ACT_annual_25.amc  0 GHz  1000 GHz  10 MHz  $((i)) deg  $CerroAPEXQ1 >data/CerroAPEX/25/ACT_annual_25.$((i)).out
    am data/ACT_annual_50.amc  0 GHz  1000 GHz  10 MHz  $((i)) deg  $CerroAPEXQ2 >data/CerroAPEX/50/ACT_annual_50.$((i)).out
    am data/ACT_annual_75.amc  0 GHz  1000 GHz  10 MHz  $((i)) deg  $CerroAPEXQ3 >data/CerroAPEX/75/ACT_annual_75.$((i)).out

    #CerroConfig Calculations
    am data/ACT_annual_25.amc  0 GHz  1000 GHz  10 MHz  $((i)) deg  $CerroConfigQ1 >data/CerroConfig/25/ACT_annual_25.$((i)).out
    am data/ACT_annual_50.amc  0 GHz  1000 GHz  10 MHz  $((i)) deg  $CerroConfigQ2 >data/CerroConfig/50/ACT_annual_50.$((i)).out
    am data/ACT_annual_75.amc  0 GHz  1000 GHz  10 MHz  $((i)) deg  $CerroConfigQ3 >data/CerroConfig/75/ACT_annual_75.$((i)).out
done