import matplotlib.pyplot as plt
import numpy as np


def plotPercentile(percentile):
    angle = 0
    file = open("data/ACT_annual_" + str(percentile) +
                "." + str(angle) + ".out", "r")
    data = file.readlines()

    x = [float(i.split(" ")[0]) for i in data]
    y = [float(i.split(" ")[1]) for i in data]

    file.close()

    # (str(percentile) + "th percentile")
    plt.plot(x, y, linewidth=1, label="Configuration")


def plotApproximate(percentile):
    angle = 0
    file = open("data/ACT_annual_" + str(percentile) +
                "_approximation." + str(angle) + ".out", "r")
    data = file.readlines()

    x = [float(i.split(" ")[0]) for i in data]
    y = [float(i.split(" ")[1]) for i in data]

    file.close()

    plt.plot(x, y, linewidth=1, label="Approximation")


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
