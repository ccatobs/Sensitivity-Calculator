import yaml
import numpy as np

a = np.array([1, 2])
print(a)


stream = open("input.yaml", 'r')
dictionary = yaml.load(stream)
for key, value in dictionary.items():
    print(key + " : " + str(value))
