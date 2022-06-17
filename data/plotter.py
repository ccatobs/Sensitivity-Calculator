import matplotlib.pyplot as plt
import numpy as np

# Values for trop_h2o_scale_factor approximations in am data generation (based on PWV ratios):
# 1st Quartile: 0.47058824 = 0.4 mm / 0.85 mm
# 3rd Quartile: 2.16470588 = 1.84 mm / 0.85 mm


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


# Plots MAM at 0 degrees Q1 using version 10.0 and 12.0
def plotVersionComparison():
    file = open("data/ACT_MAM_25.0.out", "r")
    data = file.readlines()
    x = [float(i.split(" ")[0]) for i in data]
    y = [float(i.split(" ")[2]) for i in data]
    file.close()
    plt.plot(x, y, linewidth=1, label="AM 12.0")

    file = open("data/ACT_MAM_50_pwv0.51.out", "r")
    data = file.readlines()
    x = [float(i.split(" ")[0]) for i in data]
    y = [float(i.split(" ")[2]) for i in data]
    file.close()
    plt.plot(x, y, linewidth=1, label="AM 10.0")
    plt.title("Transmission vs Frequency MAM 0 degrees Q1")


# Plots the approximation and configuration file for a given percentile and angle
def plotApproximationAccuracy(percentile, angle):
    if percentile != 25 and percentile != 75:
        print("Configuration files don't exist for that percentile")
        raise SystemExit(1)
    if angle < 15 or angle > 75:
        print("Angle out of range")
        raise SystemExit(1)
    plotFraction(angle, percentile, "Approximation", approx=True)
    plotFraction(angle, percentile, "Configuration")
    plt.title("Configuration File Approximation Accuracy " +
              str(percentile) + "th percentile, " + str(angle) + " degrees Zenith Angle")


# Plots the transmission at a given angle using configuration files
def plotPercentiles(angle):
    if angle < 15 or angle > 75:
        print("Angle out of range")
        raise SystemExit(1)
    plotFraction(angle, 25, "25th percentile/1st quartile")
    plotFraction(angle, 50, "50th percentile/2nd quartile")
    plotFraction(angle, 75, "75th percentile/3rd quartile")
    plt.title("Transmission at Various Atmospheric Conditions at " +
              str(angle) + " degrees Zenith Angle")


# Plots the transmission at a given angle using 50th percentile configuration file and an adjusted PWV
def plotApproximations(angle):
    if angle < 15 or angle > 75:
        print("Angle out of range")
        raise SystemExit(1)
    plotFraction(angle, 25, "25th percentile/1st quartile", approx=True)
    plotFraction(angle, 50, "50th percentile/2nd quartile")
    plotFraction(angle, 75, "75th percentile/3rd quartile", approx=True)
    plt.title("Transmission at Various PWV at " +
              str(angle) + " degrees Zenith Angle")


# Plots the transmission for ACT at a given angle, both the approximations and configuration files
def plotACT(angle):
    if angle < 15 or angle > 75:
        print("Angle out of range")
        raise SystemExit(1)
    plotFraction(angle, 25, "Q1 Approximation", approx=True)
    plotFraction(angle, 25, "Q1 Configuration")
    plotFraction(angle, 50, "Q2")
    plotFraction(angle, 75, "Q3 Approximation", approx=True)
    plotFraction(angle, 75, "Q3 Configuration")
    plt.title(
        "Approximations and Configuration Files at " + str(angle) + " degrees Zenith Angle")


# Plots Steve's data for a given angle
def plotSteve(angle):
    if angle < 15 or angle > 75:
        print("Angle out of range")
        raise SystemExit(1)
    plotFraction(angle, 25, "Q1", prefix="Steve/")
    plotFraction(angle, 50, "Q1", prefix="Steve/")
    plotFraction(angle, 75, "Q1", prefix="Steve/")
    plt.title("Transmission Calculated from Steve's Method at " +
              str(angle) + " degrees Zenith Angle")


# Plots the transmission for CCAT at a given angle
def plotCCAT(angle):
    if angle < 15 or angle > 75:
        print("Angle out of range")
        raise SystemExit(1)
    plotFraction(angle, 25, "Q1", prefix="CCAT/")
    plotFraction(angle, 50, "Q1", prefix="CCAT/")
    plotFraction(angle, 75, "Q1", prefix="CCAT/")
    plt.title("Transmission at CCAT Site at " +
              str(angle) + " degrees Zenith Angle")


def plotAll(angle, percentile):
    if angle < 15 or angle > 75:
        print("Angle out of range")
        raise SystemExit(1)

    plotFraction(angle, percentile, "CCAT", prefix="CCAT/")
    plotFraction(angle, percentile, "Steve's Method", prefix="Steve/")
    plotFraction(angle, percentile, "ACT")
    if angle == 45:
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
    plt.title("Transmission at Various Sites at " +
              str(angle) + " degrees Zenith Angle and Q" + str(int(percentile/25)))


plotAll(44.5, 75)
plt.ylim(ymin=0, ymax=1)
plt.xlim(xmin=0, xmax=1000)
plt.grid(which="both", axis="y")
plt.legend(loc="upper right")
plt.xlabel("Frequency (GHz)")
plt.ylabel("Transmission")
plt.show()
