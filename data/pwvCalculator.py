import numpy as np


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
    return configPWVHelper("data/" + str(percentile) + "/ACT_annual_" + str(percentile) + ".45.err")


def apexPWVHelper(filePath, percentile):
    file = open(filePath, "r")
    data = file.readlines()

    x = [float(i.split(",")[2]) if i.startswith("2") else -1 for i in data]

    file.close()
    return np.percentile(x[1:], percentile)


def apexPWV(percentile):
    return (apexPWVHelper("data/APEX PWV/2019.csv", percentile) + apexPWVHelper("data/APEX PWV/2020.csv", percentile) + apexPWVHelper("data/APEX PWV/2021.csv", percentile)) / 3


def apexYearly(percentile):
    return (apexPWVHelper("data/APEX PWV/2019.csv", percentile), apexPWVHelper("data/APEX PWV/2020.csv", percentile), apexPWVHelper("data/APEX PWV/2021.csv", percentile))


print(configPWV(25), apexPWV(25), apexYearly(25))
print(configPWV(50), apexPWV(50), apexYearly(50))
print(configPWV(75), apexPWV(75), apexYearly(75))
