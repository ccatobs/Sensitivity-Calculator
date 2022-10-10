#To use this, you must be in the Excel Spreadsheet folder and have the am code installed. It takes approximately 10 minutes to run, and lists percentage completed to stdout. The files should only be created inside folders inside data/

#! usr/bin/bash

#Values taken from https://arxiv.org/pdf/2007.04262.pdf
CerroChajnantorPWVQ1=0.36
CerroChajnantorPWVQ2=0.67
CerroChajnantorPWVQ3=1.28

#Create base files so ACT's PWVs are known
am data/ACT_annual_50.amc  0 GHz  1000 GHz  10 MHz  45 deg  1.0 >/dev/null 2>data/ACTPWV/ACT_annual_50.45.err

#PWV values in ACT configuration files
ACTConfigQ2=$(python3 pwvCalculator.py 50)

#Ratio of PWVs used in CCAT site data generation
CerroConfigQ2=$(bc <<<"scale=10; $CerroChajnantorPWVQ2/$ACTConfigQ2")


for i in {15..75}
do
    for s in {1..40}
    do
        am data/ACT_annual_50.amc  0 GHz  1000 GHz  10 MHz  $((i)) deg  $(bc <<<"scale=10; $CerroConfigQ2/20*$((s))") >data/VariablePWV/ACT_annual_$((s)).$((i)).out 2>/dev/null
    done
    PERCENT=$(bc <<<"scale=0; ($((i))-14)/0.61")
    SIGN="%"
    echo "${PERCENT}${SIGN}"
done