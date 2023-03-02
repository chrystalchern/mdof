'''Frequency Domain Decomposition (FDD), including Peak Picking (PP) technique.

This script computes the SVD from the cross PSD, obtained from the raw time series data.
 Then, a PP algorithm is used to automatically pick the modal eigenfrequencies peaks from the first singular vector,
 provided they are not too close to each other. Peaks are selected according to their magnitude, in descending steps.
'''
import pandas as pd
import numpy as np
from scipy import signal

import inputs_fdd as i_f

########## OPEN ACCELEROGRAMS FILE
df0 = []
try:
    df0 = pd.read_csv(i_f.folder_inp + i_f.file_inp + '.csv', skiprows=i_f.h_cor, header=None)
    len_df0 = len(df0)
    if len_df0 <= 1:
        raise Exception("File df0 with issues - Error: it seems empty")
except Exception as e:
    print("Problems opening event file {} - Error: {}".format(i_f.file_inp, e))

########## CALCULATE PSD
psd = []    # CPSD
#review npers
if i_f.npers > len_df0:
    i_f.npers = len_df0
    print('npers adjusted!')
nfft = i_f.npers//2+1
len_c_df0 = len(df0.columns)
psd3 = np.zeros((len_c_df0,len_c_df0,nfft)) # 3D array "(j,i,freqs_1)"

for j in range(len_c_df0):
    df0_cj = df0.iloc[:,j]
    df0_cj.reset_index(inplace=True, drop=True)

    #PSD - Cross Power Spectral Density - using Welch’s method.
    for i in range(len_c_df0):
        df0_ci = df0.iloc[:,i]
        df0_ci.reset_index(inplace=True, drop=True)

        # Frequencies and Magnitude PSD
        freqs_1, psd00 = signal.csd(df0_cj, df0_ci, fs=i_f.Freq, nperseg=i_f.npers)
        psd0 = np.abs(psd00)
        #print(freqs_1)
        if j == 0 and i == 0:
            psd.append(freqs_1)
        psd.append(psd0)
        #print(psd)
        for k in range(len(psd0)):
            psd3[j,i,k] = psd0[k]

print("Shape:", psd3.shape)

########## CALCULATE SVD
s1 = []   # 1st svd
for i in range(psd3.shape[2]):
    df3_z = psd3[:,:,i]
    U, s, Vh = np.linalg.svd(df3_z)
    s1.append(s[0])

########## FIND PEAKS
# Frequencies
Ts = 1/i_f.Freq
# Peaks limited between min and max pre-defined range
for k in range(len(freqs_1)):
    aux_f = float(freqs_1[k])
    if aux_f <= i_f.freq_r_L: #range LOW
        f_min = k
    if aux_f <= i_f.freq_r_H: #range HIGH
        f_max = k
frq_mm = freqs_1[range(f_min,f_max+1)]
s1_ = np.array(s1)
s1_mm = s1_[range(f_min,f_max+1)]

# Normalize s1
s1_norm = (s1_mm - np.min(s1_mm))/np.ptp(s1_mm)

# Find max peak
peakY = np.max(s1_norm)  #max amplitude
if np.isnan(peakY) == False:
    locY = np.where(s1_norm == peakY)[0][0]  #peak location
    freY = "%.3f" % frq_mm[locY]   #peak frequency value - three decimal places
print("Max Peak Freq:", freY)
print("position:", locY)

frqY = []
frqY.append(freY)

# Find other peaks
peaks, _ = signal.find_peaks(s1_norm, distance = i_f.dist)
aux_p = 0
harm = False
#List of all peaks - importance descending in steps according to p_peak (percentages from main peak)
p_peak = [1,0.9,0.8,0.7,0.6,0.5,0.4,0.3]
for k0 in range(1,len(p_peak)):
    for k in range(0,len(peaks)):
        frqYz = "%.3f" % frq_mm[peaks[k]]      # Get the actual frequency value
        #if it larger in the next step of importance, and if it is between range LOW and range HIGH
        if (s1_norm[peaks[k]] >= p_peak[k0]*peakY and s1_norm[peaks[k]] < p_peak[k0-1]*peakY and 
        float(frqYz) >= i_f.freq_r_L and float(frqYz) <= i_f.freq_r_H):
            #Ignore harmonics and append Frequency value to list
            harm = False
            for f0 in frqY:
                for hf in range(1,6): #not equal nor one of first 5 harmonics:
                    if float(frqYz) == float(f0)*hf:
                        harm = True
            if harm == False:
                aux_p = aux_p + 1
                if aux_p <= i_f.limit_frq-1:
                    frqY.append(frqYz)

# Save Results
df1 = pd.DataFrame(frqY, columns=['Peaks'])
print(df1)

try:
    df1.to_csv(i_f.folder_out + 'res_Peaks_' + i_f.file_inp + '.csv',index=False)
except Exception as e:
    print('Problem saving file - Error:', e)
