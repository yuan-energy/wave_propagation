import numpy as np 
import matplotlib.pyplot as plt
import sys

filename = sys.argv[1];
print("Ploting file : " + filename)
result = np.loadtxt(filename)

time = result[: ,0]
ampl = result[:,1]
plt.plot(time, ampl)
plt.show()

