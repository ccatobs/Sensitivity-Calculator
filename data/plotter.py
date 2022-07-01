import matplotlib.pyplot as plt
import numpy as np


def plotCustom(filePath, label):
    file = open(filePath, "r")
    data = file.readlines()

    x = [float(i.split(" ")[0]) for i in data]
    y = [float(i.split(" ")[1]) for i in data]

    file.close()

    plt.plot(x, y, linewidth=1, label=label)


def plotFraction(angle, percentile, label, prefix="", approx=False):
    if prefix == None:
        prefix = ""
    if angle >= 15 and angle <= 75 and int(angle) == angle:
        plotCustom("data/" + prefix + str(percentile) + ("_approx" if approx else "") + "/ACT_annual_" +
                   str(percentile) + "." + str(angle) + ".out", label)
    elif int(angle) != angle:
        floorFile = open("data/" + prefix + str(percentile) + ("_approx" if approx else "") + "/ACT_annual_" +
                         str(percentile) + "." + str(int(np.floor(angle))) + ".out", "r")
        floorData = floorFile.readlines()
        floorx = [float(i.split(" ")[0]) for i in floorData]
        floory = [float(i.split(" ")[1]) for i in floorData]
        floorFile.close()

        ceilFile = open("data/" + prefix + str(percentile) + ("_approx" if approx else "") + "/ACT_annual_" +
                        str(percentile) + "." + str(int(np.ceil(angle))) + ".out", "r")
        ceilData = ceilFile.readlines()
        ceily = [float(i.split(" ")[1]) for i in ceilData]
        ceilFile.close()

        prop = angle - np.floor(angle)
        y = [flooryi * (1 - prop) + ceilyi *
             prop for (flooryi, ceilyi) in zip(floory, ceily)]

        plt.plot(floorx, y, linewidth=1, label=label)
    else:
        print("Angle out of range")
        raise SystemExit(1)


# Plots a piecewise function for transmission values from the excel sheet
def plotExcel(percentile):
    excelTrans = [
        [.441, .253, .085],
        [.789, .679, .506],
        [.882, .813, .693],
        [.95, .923, .871],
        [.964, .946, .91],
    ]
    excelFreq = [8.5e+11, 4.05e+11, 3.48e+11, 2.85e+11, 2.22e+11]
    excelBW = [9.7e+10, 3.e+10, 3.6e+10, 7.e+10, 5.6e+10]
    for (f, t, bw) in zip(excelFreq, excelTrans, excelBW):
        x = [(f-bw/2)/1e9, (f+bw/2)/1e9]
        y = [t[int(percentile/25)-1], t[int(percentile/25)-1]]
        plt.plot(x, y,
                 linewidth=2, label=("Excel " + str(int(f/1e9)) + " GHz"))


# Plots the current transmission graphs, and excel sheet if angle is 45 degrees and excelTrans is True
def plotCurrent(angle, percentile, label=None, excelTrans=False):
    if angle < 15 or angle > 75:
        print("Angle out of range")
        raise SystemExit(1)
    if label == None:
        label = str(angle) + " degrees, Q" + str(int(percentile/25))
    plotFraction(angle, percentile, label, prefix="CerroConfig/")
    if angle == 45 and excelTrans:
        plotExcel(percentile)
    plt.title("Transmission at " +
              str(angle) + " degrees Zenith Angle, Q" + str(int(percentile/25)))


# Plots all transmission graphs
def plotAll(angle):
    plotCurrent(angle, 25, "Q1")
    plotCurrent(angle, 50, "Q2")
    plotCurrent(angle, 75, "Q3")
    plt.title("Transmission at " + str(90 - angle) +
              " Observation Elevation Angle")


# Requested graphs go here
plotAll(45)


plt.ylim(ymin=0, ymax=1)
plt.xlim(xmin=0, xmax=1000)
plt.grid(which="both", axis="y")
plt.legend(loc="upper right")
plt.xlabel("Frequency (GHz)")
plt.ylabel("Transmission")
plt.show()
