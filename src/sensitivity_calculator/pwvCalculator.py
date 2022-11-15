import sys
import os
absolute_path = os.path.dirname(__file__)


def configPWVHelper(filePath):
    file = open(filePath, "r")
    data = file.readlines()

    curr = None

    for line in data:
        if " um_pwv)" in line:
            curr = float(line.split(" um_pwv)")[0].split("(")[1])

    file.close()
    return curr/1e3


def configPWV(percentile):
    return configPWVHelper(os.path.join(absolute_path, "data/ACTPWV/ACT_annual_") + str(percentile) + ".45.err")


def datagen(percentile):
    return configPWVHelper("src/sensitivity_calculator/data/ACTPWV/ACT_annual_" + str(percentile) + ".45.err")


if __name__ == '__main__':
    print(datagen(sys.argv[1]))
