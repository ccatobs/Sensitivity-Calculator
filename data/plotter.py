import matplotlib.pyplot as plt
import numpy as np


def plotPercentile(percentile):
    angle = 45
    file = open("data/ACT_annual_" + str(percentile) +
                "." + str(angle) + ".out", "r")
    data = file.readlines()

    x = [float(i.split(" ")[0]) for i in data]
    y = [float(i.split(" ")[1]) for i in data]

    file.close()

    plt.plot(x, y, linewidth=1, label=(str(percentile) + "th percentile"))


def plotApproximate(percentile):
    angle = 45
    file = open("data/ACT_annual_" + str(percentile) +
                "_approximation." + str(angle) + ".out", "r")
    data = file.readlines()

    x = [float(i.split(" ")[0]) for i in data]
    y = [float(i.split(" ")[1]) for i in data]

    file.close()

    plt.plot(x, y, linewidth=1, label=(str(percentile) + "th percentile"))


plotPercentile(25)
# plotApproximate(25)
plotPercentile(50)
# plotApproximate(75)
plotPercentile(75)

plt.ylim(ymin=0, ymax=1)
plt.xlim(xmin=0, xmax=1000)
plt.grid(which="both", axis="y")
plt.legend(loc="upper right")
plt.title("Transmission vs Frequency 3rd Quartile")
plt.xlabel("Frequency (GHz)")
plt.ylabel("Transmission")
plt.show()
