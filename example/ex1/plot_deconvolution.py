import numpy as np 
import matplotlib.pyplot as plt
import sys

filename = sys.argv[1];
print("Ploting file : " + filename)
result = np.loadtxt(filename)

surface_wave_filename = "scaled_northridge_acc.dat"
surface_wave = np.loadtxt(surface_wave_filename)

time = result[: ,0]
ampl = result[:,1]
surface_t = surface_wave[:,0]
surface_a = surface_wave[:,1]
plt.plot(time, ampl, linewidth=2, label="Deconvolution")
plt.plot(surface_t, surface_a, linewidth=2, label="Surface")
plt.legend()
plt.xlabel("Time [second] ")
plt.ylabel("Acceleration [m/s^2] ")
plt.grid()
plt.savefig("deconvolution.jpg")
plt.show()

