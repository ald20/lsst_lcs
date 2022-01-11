import numpy as np
import pandas as pd
from scipy import stats
from astropy.timeseries import LombScargle
from gatspy.periodic import LombScargleFast
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter
import alive_progress
import os

path_67p=os.path.expanduser('~/')+'Documents/year1/shape_modelling/67p/'
writeLCsToFile=False
limiting_rh = 0 # AU, pre-defined in filename, not here.
"""
# Read in unshuffled lightcurve data
# Check which LC data you are reading in!
#lightcurveRaw = open(f'{path_67p}67P_DATA_5SIG_{limiting_rh}AU_FINAL.txt')
lightcurveRaw = open(f'{path_67p}67P_LCs_ALL_R.dat')

data = pd.read_csv(lightcurveRaw, sep=' ', header=0)
print(data.columns)
#print(data['appMag'])
mjd = data['mjd']; mag0=data['appMag']
unc=data['unc']; rh=data['rh']
delta=data['delta']; alpha=data['alpha']

print(f'Number of points in this lightcurve (>{limiting_rh} au) = {len(mag0)}')
#mag0 = np.array(mag0); unc=np.array(unc);
#alpha=np.array(alpha); delta=np.array(delta);
#rh=np.array(rh)

def correctForDistance(m, d, r):
    H_11a = m-(5.*np.log10(d*r))
    return H_11a
def correctForPhase(H_11a, b, a):
    H_110 = H_11a - (b*a)
    return H_110

# Want to randomise magnitudes from normal distribution with mean=mag and sigma=unc
# Create datatable to store randomised magnitudes in columns
# Every column is an instance of the complete randomised lightcurve
# n columns
n = 500
container = np.empty([len(mag0), n])

for i in range(len(container[1])):
    for j in range(len(mag0)):
        random_mag = np.random.normal(loc=mag0[j], scale=unc[j], size=None)
        container[j,i]=random_mag


# Now, need to perform linear regression on each column of container
# With phase angles, to guess best fitting slope (beta)
# First, need to distance correct magnitudes

abs_mags = np.empty_like(container)
for i in range(len(container[:,0])):
    # Correct every row with the same observing geometry
    for j in range(len(container[0])):
        abs_mags[i,j] = correctForDistance(container[i,j], rh[i], delta[i])

# Make array to store beta values and standard errors in
beta_array = np.empty([len(abs_mags[0]),2])

# Linear regression for each column in turn:
for k in range(len(abs_mags[0])):
    beta, h0, rval, pval, se = stats.linregress(alpha,abs_mags[:,k])
    beta_array[k]=[beta,se]
    # Use beta value to correct for phase function
    abs_mags[:,k] = correctForPhase(abs_mags[:,k], beta_array[k,0], alpha)
#print(beta_array)

float_formatter = "{:.6f}".format
np.set_printoptions(formatter={'float_kind':float_formatter})


if writeLCsToFile==True:
    # Store as a pandas dataframe for ease of printing:
    absMagsDf = pd.DataFrame(data=abs_mags, index=None, columns=None)

    lcFileName = f'{path_67p}67p_lcs_{str(n)}_{limiting_rh}au_h0.txt'
    #fout = open(lcFileName, 'w+')
    pd.set_option("display.max_rows", len(abs_mags), "display.max_columns", n)
    #print(absMagsDf)
    absMagsDf.to_csv(lcFileName, header=None, index=None, sep=' ', mode='w+', float_format='%.6f')
    #fout.close()
else:
    print('Not saving anything to file today :) ')
"""
# Plot random data phased to actual period to see if it's hot garbage or not
#p_real = 12.4041

#phased_times = (mjd-mjd[0])*24./p_real % 1

#plt.scatter(phased_times,abs_mags[:,0],c=mjd, cmap='viridis_r')
#plt.colorbar()
#plt.show()
"""
# So now we have a table of absolutely-corrected magnitudes, one lihgtcurves per column.
# Want to do a Lomb-Scargle periodogram search on each column

# period range to sample:
min_prot=0.1#*u.day # days
max_prot=2.#*u.day # days
Nsample = 1000000

period_array = np.empty([len(abs_mags[0])])
# Define frequency range:
freq = np.linspace(1./max_prot, 1./min_prot, Nsample)
periods_for_plotting = np.arange(min_prot, max_prot, 1./Nsample)
print(f'Performing LS period search on {n} randomised lightcurves:')
for i in range(len(abs_mags[0])):
    print(f'Completed {round((float(i)/n)*100.,2)}%')
    ##  Astropy periodogram implementation:
    #LSpower = LombScargle(mjd, abs_mags[:,i], unc).power(freq, method='fast')
    #best_freq = freq[np.argmax(LSpower)]
    #best_period = 2.*24.*(1./best_freq)
    ## Gatspy periodogram:
    model = LombScargleFast(optimizer_kwds={'quiet':True}).fit(mjd, abs_mags[:,i], unc)
    model.optimizer.period_range = (min_prot, max_prot)
    P,S = model.find_best_periods(n_periods=1, return_scores=True)
    #print("P=", P*24.)
    #print("S=", S)

    #print("LSpower best period=", best_period*24.)

    #scores=model.score_frequency_grid(1./max_prot, ((1./min_prot)-(1./max_prot))/Nsample, Nsample)
    #plt.plot(1./freq, scores)
    #plt.show()

    #LSpower = LombScargleFast(fit_period=True)
    #LSpower.optimizer.set(period_range=(min_prot, max_prot), first_pass_coverage=100)
    #LSpower.fit(mjd, abs_mags[:,i],unc)
    #best_period=LSpower.best_period
    #print(best_period)
    #print("Best period = {0:.4f}h".format(best_period*24.*2.))
    period_array[i] = P*24.*2.
"""
"""print(period_array)

plt.clf()
plt.hist(period_array,bins=10, range=(min(period_array), max(period_array)))

#plt.tick_params(axis='x')
plt.gca().xaxis.set_major_formatter(StrMethodFormatter('{x:,.4f}'))
plt.show()


plt.hist(beta_array[:,0], bins=10, range=(min(beta_array[:,0]), max(beta_array[:,0])))
#plt.tick_params(axis='x')
plt.gca().xaxis.set_major_formatter(StrMethodFormatter('{x:,.4f}'))
plt.show()"""



