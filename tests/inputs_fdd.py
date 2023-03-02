#Input data to peak_picking.py

folder_inp = '.'            # input folder name
file_inp = 'myshakedata.csv'         # data input file name
folder_out = 'fddresults/'             # output folder name
h_cor = 0                   # Do not consider first h_cor rows (Header)
Freq = 50                   # Sample Frequency [Hz]
freq_r_L = 0.2              # Frequency range Low value [Hz]
freq_r_H = 20               # Frequency range High value [Hz]
dist = 4                    # Min distance between peaks - samples between neighbouring peaks (dist * Freq / n = interval_Hz)
limit_frq = 10              # Maximum amount of frequencies to present
npers = 1024*2              # Length of each segment
