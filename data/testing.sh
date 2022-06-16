#To use this, you MUST be in the Excel Spreadsheet folder and have the am code (I used version 12.0) installed
#This script takes about 15 minutes on my machine; it produces console outputs via am code's automatic outputs as a sign it's running
#! usr/bin/bash


#Values taken from https://arxiv.org/pdf/2007.04262.pdf, feel free to update
CerroChajnantorPWVQ1="0.36"
CerroChajnantorPWVQ2="0.67"
CerroChajnantorPWVQ3="1.28"

#Values approximated by taking Chajnantor Plateau values from https://arxiv.org/pdf/2007.04262.pdf
ActPWVQ1="0.60"
ActPWVQ2="1.05"
ActPWVQ3="1.98"

CCATQ1=$(bc <<<"scale=2; $CerroChajnantorPWVQ1/$ActPWVQ1")
CCATQ2=$(bc <<<"scale=2; $CerroChajnantorPWVQ2/$ActPWVQ2")
CCATQ3=$(bc <<<"scale=2; $CerroChajnantorPWVQ3/$ActPWVQ3")

echo $CCATQ1
echo $CCATQ2
echo $CCATQ3