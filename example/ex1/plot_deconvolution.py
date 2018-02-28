import numpy as np 
import matplotlib.pyplot as plt
import sys

# ===========================================
# Font size
# ===========================================
import matplotlib
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 22}
matplotlib.rc('font', **font)
# ===========================================

filename = sys.argv[1];
print("Ploting file : " + filename)
result = np.loadtxt(filename)

surface_wave_filename = "scaled_northridge_acc.dat"
surface_wave = np.loadtxt(surface_wave_filename)

time = result[: ,0]
ampl = result[:,1]
surface_t = surface_wave[:,0]
surface_a = surface_wave[:,1]

plt.figure(figsize=(16,12))
plt.plot(time, ampl, 'r', linewidth=2, label="Deconvolution to Bedrock")
plt.plot(surface_t, surface_a, 'k', linewidth=2, label="Surface")
plt.legend()
plt.xlabel("Time [second] ")
plt.ylabel("Acceleration [m/s^2] ")
plt.grid()
plt.savefig("deconvolution.pdf", transparent=True, bbox_inches='tight')
plt.show()

