# CCAT-Prime Sensitivity Calculator

## Warning:

This is the development branch. For the lastest guarenteed working version, go to the "basic" branch.

## Usage:

Clone repo to local machine. Input parameters into input.yaml. Run Excel Spreadsheet.py. Output is given in output.yaml. Requires: Python 3.9.12, numpy 1.21.5, texttable 1.6.4, matplotlib 3.5.2.
If data needs to be regenrated (such as due to updated PWV measurements), then run data/amDataGenerator.sh (with updated PWV values) with GNU bc 1.07.1 and am 12.0

## Transmission Graphing:

To better visualize the calculated transmissions, run data/plotter.py. The file includes many graphs; to see each one, call the requested function near the bottom of the file around the comment "Requested graphs go here".

## Methods:

### Atmospheric Transmission:

Calculating atmospheric transmission makes use of a look-up table of transmission values calculated using AM code 12.0. AM gives configuration files for the ACT site, but not CCAT. However, the files have a parameter for scaling tropospheric water vapor. As per pages 27-28 of the AM 12.0 manual, adjusting the scale factor can approximate other site conditions. Thus, we use the ratio of PWV found in Cerro Chajnantor to the PWV in the configuration values as the tropospheric water vapor scaling parameter. Transmissions are calculated for integer degree angles between 15 and 75 degree observation elevation angles, and then interpolated to allow a continuous range of observation elevation angle inputs.

### Spill Efficiency

Spill efficiency is approximated by taking a beam profile for specific frequency and pixel spacing information and scaling the x-axis to approximate other frequency and pixel spacing information. The amount the x-axis is scaled is the ratio between the frequencies times the ratio between the pixel spacings. To use a new beam profile, change getSpillEfficiency in Excel Spreadsheet.py to load the new beam's data instead.

## Changes from Excel Sheet:

Power per pixel uses t(cold) from the corresponding frequencies instead of 850 GHz.
EoR t(cold) uses t(cold) from the corresponding frequencies instead of 405 GHz.
EoR Accepted Modes uses spill efficiency from corresponding frequencies instead of 1071um using 861um's spill efficiency
