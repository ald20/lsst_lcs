# lsst_lcs

Scripts used for analysis in LSST comet lightcurves project. Brief summary of each script:

i) Format_LSST_LCC.ipynb:
- reads columns from simulator output, parses Horizons for Earth x,y,z vectors (only sun provided by LSST). 
- Performs light time correction using Earth distance from LSST simulator output.
- Also contains MikkoWrite function, writing JDs, intensities (dummy) and vectors to Mikko file format.

ii) mag_from_phase_func.py
- Reads in phase angles from 67P_PHASE_ALL_lsst.txt (column of phase angles from LSST output)
- Calculates H(1,1,alpha) using beta and H0 from Lowry 2012
- Reads in the Matlab computed shifts from 67P_matlab_LC_shifts.txt
- ADDS the shift to its corresponding H(1,1,alpha)
- Writes the shifted magnitude to 67P_mag_list_ALL.txt
- Also generates intensities for file 67P_int_list_ALL.txt (using corr_fac as the mean of the shifted magnitudes)

iii) LSST_LC_procedures.ipynb:
- Reads in magnitudes for use (these need to be created in first place) from 67p_20230109_I11_R.dat (MJD, absMag, SunXyz, r_h, EarthXyz, Delta)
- Calculates an apparent magnitude m: H = m - 5log(delta*rh)
- **Uncertainties prepared here** function created from manual FORS2 SNR inputs
- Writes columns mjd, mag, m_app, rh, delta, alpha, snr, unc to 67P_LCs_ALL_R.dat
