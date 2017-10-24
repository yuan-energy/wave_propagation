import numpy as np 
import matplotlib.pyplot as plt
import sys

# filename = sys.argv[1];
# print("Ploting file : " + filename)
# result = np.loadtxt(filename)

data = []
step_len = 25
for x in range(0,251,step_len):
	wave_filename = "wave_at_depth_"+str(x)+"_acc.txt"
	wave = np.loadtxt(wave_filename)
	data.append(wave[:,1])

wave_filename = "ricker_acc.dat"
wave = np.loadtxt(wave_filename)
time = wave[: ,0]

for x in range(len(data)) :
	layer = data[x]
	layer = [item + (24 -x)*1.6 for item in layer ]
	plt.plot(time, layer, linewidth=2, label= "depth" + str(x*step_len) )

plt.legend()
plt.grid()
plt.yticks([])
plt.xlabel("Times/(s)")
plt.savefig("upward_wave.jpg")
plt.show()




