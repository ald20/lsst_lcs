# lsst_lcs

Scripts used for analysis in LSST comet lightcurves project. Brief summary of each script:

i) Format_LSST_LCC.ipynb:
- reads columns from simulator output (LSST_out_LC_67P.txt), parses Horizons for Earth x,y,z vectors (only sun provided by LSST). 
- Performs light time correction using Earth distance from LSST simulator output.
- Also contains MikkoWrite function, writing JDs, intensities (dummy) and vectors to Mikko file format.
- Gets solar phase angles and ecliptic longitudes, latitudes from Horizons ephemerides.
- Final steps in the script: reads in shifted absolute mags (from py script below) and makes MikkoFile with real (non-dummy) intensities (67P_LC_NU.dat)
- Also makes 20230109_I11_R.dat file with unc 0.01 mag in third column

ii) mag_from_phase_func.py
- Reads in phase angles from 67P_PHASE_ALL_lsst.txt (column of phase angles from LSST output)
- Calculates H(1,1,alpha) using beta and H0 from Lowry 2012
- Reads in the Matlab computed shifts from 67P_matlab_LC_shifts.txt
- ADDS the shift to its corresponding H(1,1,alpha)
- Writes the shifted magnitude to 67P_mag_list_ALL.txt
- Also generates intensities for file 67P_int_list_ALL.txt (using corr_fac as the mean of the shifted magnitudes)

iii) LSST_LC_procedures.ipynb:
- Reads in magnitudes for use (these need to be created in first place) from 67p_20230109_I11_R.dat (MJD, absMag, SunXyz, r_h, EarthXyz, Delta, EclLon, EclLat)
- Calculates an apparent magnitude m: H = m - 5log(delta*rh)
- **Uncertainties prepared here** function created from manual FORS2 SNR inputs (saved in file FORS2_snr_vals.dat)
- Writes columns mjd, mag, m_app, rh, delta, alpha, snr, unc to 67P_LCs_R.dat

iv) periodogram.ipynb
- Reads in data contained in table 67P_LCs_R.dat
- Treating the apparent magnitudes here as though they are real, calibrated data points, convert them to H(1,1,alpha) using an 'average' phase function.
- Parses Horizons to find next perihelion MJDs for 67P
- Performs periodogram analysis on data points (option to limit by heliocentric distance), searching rotation periods between 0.1 and 2. days.

v) LSST_baselinev13_query.py
- Script to search baseline_v1.3_10years.db using SQLITE to return 5 sigma depths and observationID
