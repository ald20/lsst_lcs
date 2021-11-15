import matplotlib.pyplot as plt
import numpy as np

mag_file = open('../67p/67P_mag_list_ALL.txt', 'r')
mjd_file = open('../67p/67P_LCS_ALL_R.dat', 'r')
lines = mag_file.readlines()
mags = []
for line in lines:
    data = line.split()
    mags.append(float(data[0]))
mag_file.close()

lines = mjd_file.readlines()
mjds = []
for line in lines[1:]:
    data = line.split()
    mjds.append(float(data[0]))

alphas=[]
alpha_file = open('../67p/67P_PHASE_ALL_lsst.txt', 'r')
lines = alpha_file.readlines()
for line in lines:
    data = line.split()
    alphas.append(float(data[0]))
alphas = np.array(alphas)
alpha_file.close()

mjds=np.array(mjds)
period=12.4041
h0 = 15.43
beta = 0.059
phase = mjds / (period*24.) % 1
print(np.mean(mags))

plt.scatter(phase,(mags+(beta*alphas)),c=mjds,cmap='viridis_r')
plt.gca().invert_yaxis()
plt.show()
