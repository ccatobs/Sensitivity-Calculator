import matplotlib.pyplot as plt
import numpy as np

# Values for trop_h2o_scale_factor approximations in am data generation (based on PWV ratios):
# 1st Quartile: 0.47058824 = 0.4 mm / 0.85 mm
# 3rd Quartile: 2.16470588 = 1.84 mm / 0.85 mm


def plotPercentile(percentile, angle, label, season):
    filePath = None
    if season == None:
        filePath = "data/" + str(percentile) + "/ACT_annual_" + \
            str(percentile) + "." + str(angle) + ".out"
    else:
        filePath = "data/" + season + "/" + \
            str(percentile) + "/ACT_" + season + "_" + \
            str(percentile) + "." + str(angle) + ".out"
    file = open(filePath, "r")
    data = file.readlines()

    x = [float(i.split(" ")[0]) for i in data]
    y = [float(i.split(" ")[1]) for i in data]

    file.close()

    plt.plot(x, y, linewidth=1, label=label)


def plotApproximate(percentile, angle, label, season):
    filePath = None
    if season == None:
        if percentile == 50:
            filePath = "data/50/ACT_annual_" + \
                str(percentile) + "_approximation." + str(angle) + ".out"
        else:
            filePath = "data/" + str(percentile) + "_approx/ACT_annual_" + \
                str(percentile) + "_approximation." + str(angle) + ".out"
    else:
        if percentile == 50:
            filePath = "data/" + season + "/50/ACT_" + \
                season + "_50." + str(angle) + ".out"
        else:
            filePath = "data/" + season + "/" + str(percentile) + "_approx/ACT_" + season + "_" + str(
                percentile) + "_approximation." + str(angle) + ".out"
    file = open(filePath, "r")
    data = file.readlines()

    x = [float(i.split(" ")[0]) for i in data]
    y = [float(i.split(" ")[1]) for i in data]

    file.close()

    plt.plot(x, y, linewidth=1, label=label)


def parseSeason(season):
    if season == None:
        return "Annually"
    else:
        return "during season " + season

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
def plotApproximationAccuracy(percentile, angle, season):
    if percentile != 25 and percentile != 75:
        print("Configuration files don't exist for that percentile")
        raise SystemExit(1)
    if int(angle) != angle or angle < 0 or angle >= 90:
        print("Interpolation not yet implemented")
        raise SystemExit(1)
    if season != None and season != "DJF" and season != "JFM" and season != "JJA" and season != "MAM" and season != "SON":
        print("Season doesn't exist")
        raise SystemExit(1)
    plotApproximate(percentile, angle, "Approximation", season)
    plotPercentile(percentile, angle, "Configuration", season)
    plt.title("Configuration File Approximation Accuracy " +
              str(percentile) + "th percentile, " + str(angle) + " degrees Zenith Angle")


# Plots the transmission at a given angle using configuration files
def plotPercentiles(angle, season):
    if int(angle) != angle or angle < 0 or angle >= 90:
        print("Interpolation not yet implemented")
        raise SystemExit(1)
    if season != None and season != "DJF" and season != "JFM" and season != "JJA" and season != "MAM" and season != "SON":
        print("Season doesn't exist")
        raise SystemExit(1)
    plotPercentile(25, angle, "25th percentile/1st quartile", season)
    plotPercentile(50, angle, "50th percentile/2nd quartile", season)
    plotPercentile(75, angle, "75th percentile/3rd quartile", season)
    plt.title("Transmission at Various Atmospheric Conditions at " +
              str(angle) + " degrees Zenith Angle")


# Plots the transmission at a given angle using 50th percentile configuration file and an adjusted PWV
def plotApproximations(angle, season):
    if int(angle) != angle or angle < 0 or angle >= 90:
        print("Interpolation not yet implemented")
        raise SystemExit(1)
    if season != None and season != "DJF" and season != "JFM" and season != "JJA" and season != "MAM" and season != "SON":
        print("Season doesn't exist")
        raise SystemExit(1)
    plotApproximate(25, angle, "25th percentile/1st quartile", season)
    plotPercentile(50, angle, "50th percentile/2nd quartile", season)
    plotApproximate(75, angle, "75th percentile/3rd quartile", season)
    plt.title("Transmission at Various PWV at " +
              str(angle) + " degrees Zenith Angle")


# Plots the transmission at a given angle, both the approximations and configuration files
def plotAll(angle, season):
    if int(angle) != angle or angle < 0 or angle >= 90:
        print("Interpolation not yet implemented")
        raise SystemExit(1)
    if season != None and season != "DJF" and season != "JFM" and season != "JJA" and season != "MAM" and season != "SON":
        print("Season doesn't exist")
        raise SystemExit(1)
    plotApproximate(25, angle, "Q1 Approximation", season)
    plotPercentile(25, angle, "Q1 Configuration", season)
    plotPercentile(50, angle, "Q2", season)
    plotApproximate(75, angle, "Q3 Approximation", season)
    plotPercentile(75, angle, "Q3 Configuration", season)
    plt.title(
        "Approximations and Configuration Files at " + str(angle) + " degrees Zenith Angle")


plotAll(45, None)
plt.ylim(ymin=0, ymax=1)
plt.xlim(xmin=0, xmax=1000)
plt.grid(which="both", axis="y")
plt.legend(loc="upper right")
plt.xlabel("Frequency (GHz)")
plt.ylabel("Transmission")
plt.show()
