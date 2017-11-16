import numpy as np 
import matplotlib.pyplot as plt
import sys

filename = sys.argv[1];
# print("Ploting file : " + filename)
# result = np.loadtxt(filename)


wave = np.loadtxt(filename)
time = wave[: ,0]
ampl = wave[:, 1]

plt.plot(time, ampl, linewidth=2, label= filename )

plt.legend()
plt.grid()
plt.savefig(filename+".jpg")
plt.show()

