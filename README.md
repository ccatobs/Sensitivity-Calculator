# CCAT-Prime Sensitivity Calculator

## Usage:

Clone repo to local machine. Input requested parameters into input.yaml. Run Excel Spreadsheet.py. Output is given in output.yaml. Requires: Python 3.9.12, numpy 1.21.5.
If data needs to be regenrated (such as due to updated PWV measurements), then run data/amDataGenerator.sh (with updated PWV values) with GNU bc 1.07.1 and am 12.0

## Transmission Graphing:

To better visualize the calculated transmissions, run data/plotter.py. The file includes many graphs; to see each one, call the requested function near the bottom of the file around the comment "Requested graphs go here". Requires: Python 3.9.12, matplotlib 3.5.2

## Methods:

The currently used method for calculating atmospheric transmission makes use of a look-up table of transmission values calculated using AM code 12.0. AM gives configuration files for the ACT site, but not CCAT. However, the files have a parameter for scaling tropospheric water vapor. As per pages 27-28 of the AM 12.0 manual, adjusting the scale factor can approximate other site conditions. Thus, we use the ratio of PWV found in Cerro Chajnantor to Chajnantor Plateau (which is at the same elevation as ACT's site) as the scale factor while using ACT's configuration files.
