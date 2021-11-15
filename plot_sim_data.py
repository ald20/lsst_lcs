import numpy as np
import matplotlib.pyplot as plt
import os
import scipy
from scipy import stats

home_dir = os.path.expanduser('~/')
path_67p = home_dir+'Documents/year1/shape_modelling/67p/'
pics_dir = path_67p+'paper_pics/'

# Read in shifted magnitudes:
shifted_mags = []
magfile = open(path_67p+'67p_mag_list_ALL.txt', 'r')
mag_lines = magfile.readlines()
for line in mag_lines:
    data = line.split()
    shifted_mags.append(float(data[0]))

# Read in phase angles
phase_angles = []
phasefile = open(path_67p+'67P_PHASE_ALL_lsst.txt', 'r')
phase_lines = phasefile.readlines()
for line in phase_lines:
    data = line.split()
    phase_angles.append(float(data[0]))

# Read in MJDs
mjd_list = []
mjdfile = open(path_67p+'67P_LCS_ALL_R.dat')
mjdLines = mjdfile.readlines()
for line in mjdLines[1:]:
    data=line.split()
    mjd_list.append(float(data[0]))

# Straight line fit to data
slope, intercept, r_value, p_value, std_err = stats.linregress(phase_angles,shifted_mags)
print(f'Slope={slope:.3f}, intercept={intercept:.3f}')

plt.scatter(phase_angles, shifted_mags, c=mjd_list, cmap='viridis_r')
plt.plot(phase_angles, (slope*np.array(phase_angles))+intercept, color = 'gray', linestyle='--', label='%.4f'%slope)
plt.legend()
plt.gca().invert_yaxis()
cbar = plt.colorbar()
plt.show()

