import matplotlib.pyplot as plt
import numpy as np

angle = 45
percentile = 25
file = open("ACT_annual_" + str(percentile) + "." + str(angle) + ".out", "r")
data = file.readlines()

x = [float(i.split(" ")[0]) for i in data]
y = [float(i.split(" ")[1]) for i in data]

file.close()

plt.plot(x, y, linewidth=1)


angle = 45
percentile = 50
file = open("ACT_annual_" + str(percentile) + "." + str(angle) + ".out", "r")
data = file.readlines()

x = [float(i.split(" ")[0]) for i in data]
y = [float(i.split(" ")[1]) for i in data]

file.close()

plt.plot(x, y, linewidth=1)


angle = 45
percentile = 75
file = open("ACT_annual_" + str(percentile) + "." + str(angle) + ".out", "r")
data = file.readlines()

x = [float(i.split(" ")[0]) for i in data]
y = [float(i.split(" ")[1]) for i in data]

file.close()

plt.plot(x, y, linewidth=1)
plt.ylim(ymin=0, ymax=1)
plt.xlim(xmin=0, xmax=1000)
plt.grid(which="both", axis="y")
plt.title("Transmission vs Frequency")
plt.xlabel("Frequency (GHz)")
plt.ylabel("Transmission")
plt.show()
