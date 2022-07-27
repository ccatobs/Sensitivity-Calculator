import sys


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


if __name__ == '__main__':
    print(configPWV(sys.argv[1]))
