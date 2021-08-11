import numpy as np
import os

home_path = os.path.expanduser('~/')
path_67p = home_path+'Documents/year1/shape_modelling/67p/'

input_file = path_67p+'67P_PHASE_ALL_lsst.txt' # List of phase angles

phase_angles = []
f = open(input_file, 'r')
lines = f.readlines()
for line in lines:
    data = line.split()
    phase_angles.append(float(data[0]))

# Define linear phase function
def phase_funcl(a):
    beta = 0.059 # mag/deg, from Lowry 2012 (\pm 0.006)
    H_R0 = 15.43
    y = (beta*a)+H_R0
    return y

# Generate distance-corrected magnitudes
mags_zero = []
for ang in phase_angles:
    mag = phase_funcl(ang)
    #mag = round(phase_funcl(ang), 4)
    mags_zero.append(mag)
    #print(mag)

# read in magnitude shifts from Matlab
shift_list = []
shifts_file = path_67p+'67P_matlab_LC_shifts.txt'
sfile = open(shifts_file, 'r')
lines = sfile.readlines()
for line in lines:
    data=line.split()
    shift_list.append(float(data[0]))

# Calcu mags and save them to list
fout=open(path_67p+'67P_mag_list_ALL.txt', 'w+')
mags_shifted = []

for i,val in enumerate(shift_list):
    shifted_mag = val+mags_zero[i]
    print('%.6f %.6f %.6f'%(mags_zero[i], val, shifted_mag))
    fout.write('%.6f\n'%shifted_mag)
    mags_shifted.append(shifted_mag)
fout.close()


# Convert these magnitudes to intensities for Mikko lc
# Normalise: subtract from each value integer mean magzero
corr_fac = round(np.mean(mags_shifted),0)
#print(corr_fac)
# Make another file containing these intensities
fout = open(path_67p+'67P_int_list_ALL.txt', 'w+')
for mag in mags_shifted:
    intens = 10**((mag-corr_fac)/(-2.5))
    fout.write('%.6f\n'%intens) #should be somewhere around 1

fout.close()
print(corr_fac)