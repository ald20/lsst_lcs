import numpy as np
import os
import pandas as pd

home_path = os.path.expanduser('~/')
path_67p = home_path+'Documents/year1/shape_modelling/67p/'

input_file = path_67p+'67P_PHASE_ALL_lsst.txt' # List of phase angles

dir_str = 'period_param/' # dir in which updated param files are stored
period_str = '11.1' # value of period to extract from/append to file names

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
# Make sure you're using the right shift file!
shift_list = []
shifts_file = f'{path_67p}{dir_str}67p_{period_str}_shifts_ALL.txt'
if not os.path.exists(shifts_file):
    print(f'Shifts file not found. Check string params for file name. Currently using {period_str} h.')
    exit()
else:
    print(f'Reading shifts from {shifts_file}')

data=pd.read_csv(shifts_file, header=None)
data.columns = ["shifts"]

# Calc mags and save them to list
fout=open(f'{path_67p}{dir_str}67P_mag_list_{period_str}_ALL.txt', 'w+')
mags_shifted = []

for i,val in enumerate(data["shifts"]):
    shifted_mag = val+mags_zero[i]
    print(f'{mags_zero[i]:.6f} {val:.6f} {shifted_mag:.6f}')
    fout.write('%.6f\n'%shifted_mag)
    mags_shifted.append(shifted_mag)
fout.close()

# Convert these magnitudes to intensities for Mikko lc
# Normalise: subtract from each value integer mean magzero
corr_fac = round(np.mean(mags_shifted),0)
#print(corr_fac)
# Make another file containing these intensities
fout = open(f'{path_67p}{dir_str}67p_int_list_{period_str}_ALL.txt', 'w+')
for mag in mags_shifted:
    intens = 10**((mag-corr_fac)/(-2.5))
    fout.write('%.6f\n'%intens) #should be somewhere around 1

fout.close()
print(corr_fac)