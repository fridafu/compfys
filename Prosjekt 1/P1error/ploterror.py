from numpy import *
from matplotlib import pyplot as plt

myfile = open("error100.txt")
error = zeros(100)
i=0
for line in myfile:
    error[i] = float(line)
    i += 1

x = linspace(0,1,100)
plt.plot(x,abs(error))
plt.show()
