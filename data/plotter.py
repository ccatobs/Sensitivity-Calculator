import matplotlib.pyplot as plt
import numpy as np

# Values for trop_h2o_scale_factor approximations in am data generation (based on PWV ratios):
# 1st Quartile: 0.47058824
# 3rd Quartile: 2.16470588

approcimationAccuracyMode = False


def plotPercentile(percentile):
    angle = 0
    file = open("data/ACT_annual_" + str(percentile) +
                "." + str(angle) + ".out", "r")
    data = file.readlines()

    x = [float(i.split(" ")[0]) for i in data]
    y = [float(i.split(" ")[1]) for i in data]

    file.close()

    label = str(percentile) + "th percentile"
    if approcimationAccuracyMode:
        label = "Configuration"
    plt.plot(x, y, linewidth=1, label=label)


def plotApproximate(percentile):
    angle = 0
    file = open("data/ACT_annual_" + str(percentile) +
                "_approximation." + str(angle) + ".out", "r")
    data = file.readlines()

    x = [float(i.split(" ")[0]) for i in data]
    y = [float(i.split(" ")[1]) for i in data]

    file.close()

    label = str(percentile) + "th percentile"
    if approcimationAccuracyMode:
        label = "Approximation"
    plt.plot(x, y, linewidth=1, label=label)


# Comment/Uncomment based on intended graph
plotApproximate(25)
plotPercentile(25)
# plotPercentile(50)
# plotApproximate(75)
# plotPercentile(75)

plt.ylim(ymin=0, ymax=1)
plt.xlim(xmin=0, xmax=1000)
plt.grid(which="both", axis="y")
plt.legend(loc="upper right")
plt.title("Transmission vs Frequency 1st Quartile")
plt.xlabel("Frequency (GHz)")
plt.ylabel("Transmission")
plt.show()
