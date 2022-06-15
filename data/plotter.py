import matplotlib.pyplot as plt
import numpy as np

# Values for trop_h2o_scale_factor approximations in am data generation (based on PWV ratios):
# 1st Quartile: 0.47058824 = 0.4 mm / 0.85 mm
# 3rd Quartile: 2.16470588 = 1.84 mm / 0.85 mm


def plotPercentile(percentile, angle, label):
    file = open("data/" + str(percentile) + "/ACT_annual_" + str(percentile) +
                "." + str(angle) + ".out", "r")
    data = file.readlines()

    x = [float(i.split(" ")[0]) for i in data]
    y = [float(i.split(" ")[1]) for i in data]

    file.close()

    plt.plot(x, y, linewidth=1, label=label)


def plotApproximate(percentile, angle, label):
    filePath = "data/" + str(percentile) + "_approx/ACT_annual_" + \
        str(percentile) + "_approximation." + str(angle) + ".out"
    if percentile == 50:
        filePath = "data/50/ACT_annual_" + \
            str(percentile) + "_approximation." + str(angle) + ".out"
    file = open(filePath, "r")
    data = file.readlines()

    x = [float(i.split(" ")[0]) for i in data]
    y = [float(i.split(" ")[1]) for i in data]

    file.close()

    plt.plot(x, y, linewidth=1, label=label)


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
    if angle != 0 and angle != 45:
        print("Approximations are only computed for 0 and 45 degrees")
        raise SystemExit(1)
    plotApproximate(percentile, angle, "Approximation")
    plotPercentile(percentile, angle, "Configuration")
    plt.title("Configuration File Approximation Accuracy " +
              str(percentile) + "th percentile, " + str(angle) + " degrees Zenith Angle")


# Plots the transmission at a given angle using configuration files
def plotPercentiles(angle):
    if int(angle) != angle or angle < 0 or angle >= 90:
        print("Interpolation not yet implemented")
        raise SystemExit(1)
    plotPercentile(25, angle, "25th percentile/1st quartile")
    plotPercentile(50, angle, "50th percentile/2nd quartile")
    plotPercentile(75, angle, "75th percentile/3rd quartile")
    plt.title("Transmission at Various Atmospheric Conditions at " +
              str(angle) + " degrees Zenith Angle")


# Plots the transmission at a given angle using 50th percentile configuration file and an adjusted PWV
def plotApproximations(angle):
    if angle != 0 and angle != 45:
        print("Approximations are only computed for 0 and 45 degrees")
        raise SystemExit(1)
    plotApproximate(25, angle, "25th percentile/1st quartile")
    plotPercentile(50, angle, "50th percentile/2nd quartile")
    plotApproximate(75, angle, "75th percentile/3rd quartile")
    plt.title("Transmission at Various PWV at " +
              str(angle) + " degrees Zenith Angle")


# Plots the transmission at a given angle, both the approximations and configuration files
def plotAll(angle):
    if angle != 0 and angle != 45:
        print("Approximations are only computed for 0 and 45 degrees")
        raise SystemExit(1)
    plotApproximate(25, angle, "Q1 Approximation")
    plotPercentile(25, angle, "Q1 Configuration")
    plotPercentile(50, angle, "Q2")
    plotApproximate(75, angle, "Q3 Approximation")
    plotPercentile(75, angle, "Q3 Configuration")
    plt.title(
        "Transmission vs Frequency Approximations and Configuration Files at " + str(angle) + " degrees Zenith Angle")


plotPercentiles(80)
plt.ylim(ymin=0, ymax=1)
plt.xlim(xmin=0, xmax=1000)
plt.grid(which="both", axis="y")
plt.legend(loc="upper right")
plt.xlabel("Frequency (GHz)")
plt.ylabel("Transmission")
plt.show()
