import matplotlib.pyplot as plt
import numpy as np

angle = 45
percentile = 75
file = open("ACT_annual_" + str(percentile) + "." + str(angle) + ".out", "r")
data = file.readlines()

x = [float(i.split(" ")[0]) for i in data]
y = [float(i.split(" ")[1]) for i in data]

file.close()

plt.plot(x, y, linewidth=1)
plt.title("Transmission vs Frequency")
plt.xlabel("Frequency (GHz)")
plt.ylabel("Transmission")
plt.show()
