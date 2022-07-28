#To use this, you must be in the Excel Spreadsheet folder and have the am code installed. It takes approximately 2 minutes to run, and produces no output to stdout. The only things created should be files inside folders in data/
#! usr/bin/bash

#Values taken from https://arxiv.org/pdf/2007.04262.pdf
CerroChajnantorPWVQ1=0.36
CerroChajnantorPWVQ2=0.67
CerroChajnantorPWVQ3=1.28

#Create base files so ACT's PWVs are known
for p in 25 50 75
do
    am data/ACT_annual_$((p)).amc  0 GHz  1000 GHz  10 MHz  45 deg  1.0 >data/$((p))/ACT_annual_$((p)).45.out 2>data/$((p))/ACT_annual_$((p)).45.err
done

#PWV values in ACT configuration files
ACTConfigQ1=$(python3 pwvCalculator.py 25)
ACTConfigQ2=$(python3 pwvCalculator.py 50)
ACTConfigQ3=$(python3 pwvCalculator.py 75)

#Ratio of PWVs used in CCAT site data generation
CerroConfigQ1=$(bc <<<"scale=10; $CerroChajnantorPWVQ1/$ACTConfigQ1")
CerroConfigQ2=$(bc <<<"scale=10; $CerroChajnantorPWVQ2/$ACTConfigQ2")
CerroConfigQ3=$(bc <<<"scale=10; $CerroChajnantorPWVQ3/$ACTConfigQ3")

#Ratio using PWV of 0.01 lower
LowerQ1=$(bc <<<"scale=10; ($CerroChajnantorPWVQ1 - 0.01)/$ACTConfigQ1")
LowerQ2=$(bc <<<"scale=10; ($CerroChajnantorPWVQ2 - 0.01)/$ACTConfigQ2")
LowerQ3=$(bc <<<"scale=10; ($CerroChajnantorPWVQ3 - 0.01)/$ACTConfigQ3")

#Ratio using PWV of 0.01 higher
HigherQ1=$(bc <<<"scale=10; ($CerroChajnantorPWVQ1 + 0.01)/$ACTConfigQ1")
HigherQ2=$(bc <<<"scale=10; ($CerroChajnantorPWVQ2 + 0.01)/$ACTConfigQ2")
HigherQ3=$(bc <<<"scale=10; ($CerroChajnantorPWVQ3 + 0.01)/$ACTConfigQ3")

#Ratio using half the PWV
HalfQ1=$(bc <<<"scale=10; ($CerroChajnantorPWVQ1 / 2)/$ACTConfigQ1")
HalfQ2=$(bc <<<"scale=10; ($CerroChajnantorPWVQ2 / 2)/$ACTConfigQ2")
HalfQ3=$(bc <<<"scale=10; ($CerroChajnantorPWVQ3 / 2)/$ACTConfigQ3")

echo $ACTConfigQ1
echo $ACTConfigQ2
echo $ACTConfigQ3

echo "----------------------------------------------"
echo $HalfQ1
echo $LowerQ1
echo $CerroConfigQ1
echo $HigherQ1

echo "----------------------------------------------"
echo $HalfQ2
echo $LowerQ2
echo $CerroConfigQ2
echo $HigherQ2

echo "----------------------------------------------"
echo $HalfQ3
echo $LowerQ3
echo $CerroConfigQ3
echo $HigherQ3

for i in {15..75}
do
    #CerroConfig Calculations
    am data/ACT_annual_25.amc  0 GHz  1000 GHz  10 MHz  $((i)) deg  $CerroConfigQ1 >data/CerroConfig/25/ACT_annual_25.$((i)).out 2>data/CerroConfig/25/ACT_annual_25.$((i)).err
    am data/ACT_annual_50.amc  0 GHz  1000 GHz  10 MHz  $((i)) deg  $CerroConfigQ2 >data/CerroConfig/50/ACT_annual_50.$((i)).out 2>data/CerroConfig/50/ACT_annual_50.$((i)).err
    am data/ACT_annual_75.amc  0 GHz  1000 GHz  10 MHz  $((i)) deg  $CerroConfigQ3 >data/CerroConfig/75/ACT_annual_75.$((i)).out 2>data/CerroConfig/75/ACT_annual_75.$((i)).err

    #Lower PWV
    am data/ACT_annual_25.amc  0 GHz  1000 GHz  10 MHz  $((i)) deg  $LowerQ1 >data/Lower/25/ACT_annual_25.$((i)).out 2>data/Lower/25/ACT_annual_25.$((i)).err
    am data/ACT_annual_50.amc  0 GHz  1000 GHz  10 MHz  $((i)) deg  $LowerQ2 >data/Lower/50/ACT_annual_50.$((i)).out 2>data/Lower/50/ACT_annual_50.$((i)).err
    am data/ACT_annual_75.amc  0 GHz  1000 GHz  10 MHz  $((i)) deg  $LowerQ3 >data/Lower/75/ACT_annual_75.$((i)).out 2>data/Lower/75/ACT_annual_75.$((i)).err

    #Higher PWV
    am data/ACT_annual_25.amc  0 GHz  1000 GHz  10 MHz  $((i)) deg  $HigherQ1 >data/Higher/25/ACT_annual_25.$((i)).out 2>data/Higher/25/ACT_annual_25.$((i)).err
    am data/ACT_annual_50.amc  0 GHz  1000 GHz  10 MHz  $((i)) deg  $HigherQ2 >data/Higher/50/ACT_annual_50.$((i)).out 2>data/Higher/50/ACT_annual_50.$((i)).err
    am data/ACT_annual_75.amc  0 GHz  1000 GHz  10 MHz  $((i)) deg  $HigherQ3 >data/Higher/75/ACT_annual_75.$((i)).out 2>data/Higher/75/ACT_annual_75.$((i)).err

    #Half PWV
    am data/ACT_annual_25.amc  0 GHz  1000 GHz  10 MHz  $((i)) deg  $HalfQ1 >data/Half/25/ACT_annual_25.$((i)).out 2>data/Half/25/ACT_annual_25.$((i)).err
    am data/ACT_annual_50.amc  0 GHz  1000 GHz  10 MHz  $((i)) deg  $HalfQ2 >data/Half/50/ACT_annual_50.$((i)).out 2>data/Half/50/ACT_annual_50.$((i)).err
    am data/ACT_annual_75.amc  0 GHz  1000 GHz  10 MHz  $((i)) deg  $HalfQ3 >data/Half/75/ACT_annual_75.$((i)).out 2>data/Half/75/ACT_annual_75.$((i)).err

    PERCENT=$(bc <<<"scale=0; ($((i))-14)/0.61")
    SIGN="%"
    echo "${PERCENT}${SIGN}"
done