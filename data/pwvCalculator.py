import numpy as np


def configPWVHelper(filePath):
    file = open(filePath, "r")
    data = file.readlines()

    pwv = 0
    currPressure = None

    for line in data:
        if "Pbase " in line:
            currPressure = float(line.split("Pbase ")[1].split(" ")[0])
        if "column h2o vmr " in line:
            pwv += currPressure * float(line.split("column h2o vmr ")[1])
            currPressure = None

    file.close()
    return pwv


def configPWV(percentile):
    return configPWVHelper("data/ACT_annual_" + str(percentile) + ".amc")


def apexPWVHelper(filePath, percentile):
    file = open(filePath, "r")
    data = file.readlines()

    x = [float(i.split(",")[2]) if i.startswith("2") else -1 for i in data]

    file.close()
    return np.percentile(x[1:], percentile)


def apexPWV(percentile):
    return (apexPWVHelper("data/APEX PWV/2019.csv", percentile) + apexPWVHelper("data/APEX PWV/2020.csv", percentile) + apexPWVHelper("data/APEX PWV/2021.csv", percentile)) / 3


def custom(percentile):
    return (apexPWVHelper("data/APEX PWV/2019.csv", percentile), apexPWVHelper("data/APEX PWV/2020.csv", percentile), apexPWVHelper("data/APEX PWV/2021.csv", percentile))


print(configPWV(25), apexPWV(25), custom(25))
print(configPWV(50), apexPWV(50), custom(50))
print(configPWV(75), apexPWV(75), custom(75))
