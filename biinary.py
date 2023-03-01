import numpy as np
from pandas import DataFrame as df

n_events = 1
n_combo = 10
n_channels = 8
n_time = 20


for e in range(n_events):
    data = np.arange(n_channels*n_time).reshape(n_time, n_channels).T

    for i in range(1,n_combo):
        bin_i = bin(i)
        str_i = str(bin_i)[2:]
        lst_i = list(map(int,str_i))
        idx_i = np.array([0]*(n_channels-len(lst_i)) + lst_i, dtype=bool)

        # print(bin_i)
        # print(str_i)
        # print(lst_i)
        print(idx_i)

