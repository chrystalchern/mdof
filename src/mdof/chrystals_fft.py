import numpy as np

class TimeSeries:
    def __init__(self, series: np.ndarray, step: float) -> None:
        if not isinstance(series, np.ndarray):
            raise TypeError("series must be a 1D array.")
        if len(series.shape)>1:
            raise Exception("received an ND array. series must be a 1D array.")
        self.N = len(series)
        self.N = self.N - self.N%2
        self.series = series[:self.N]
        
        if not isinstance(step, float):
            raise TypeError("step must be a float.")
        self.dt = step

    def chrystals_dft(self):
        N = self.N    
        strangs_omega_bar = np.exp(-1j*2*np.pi/N)
        inv_fourier_matrix = strangs_omega_bar**np.outer(np.arange(N),np.arange(N))
        fourier_coeffs = 2.0/N*np.abs(inv_fourier_matrix @ self.series)[1:N//2]
        fourier_freqs = np.arange(1,N/2)/N/self.dt
        return (fourier_freqs,fourier_coeffs)

def chrystals_fft(series):
    N = len(series)
    if N == 1:
        return series
    else:
        strangs_omega_bar = np.exp(-1j*2*np.pi/N)
        fft_evens = chrystals_fft(series[0:N:2])
        fft_odds = chrystals_fft(series[1:N:2])
        D = strangs_omega_bar**np.arange(N/2)
        return np.concatenate((fft_evens+D*fft_odds, fft_evens-D*fft_odds))


if __name__ == "__main__":
    # import sys
    # chrystals_dft(sys.argv[0])

    tf = 10
    nt = 256
    time = np.linspace(start=0,stop=tf,num=nt)
    dt = time[1]-time[0]
    series = np.sin(5*2*np.pi*time)+np.random.normal(0,0.1,nt)

    from matplotlib import pyplot as plt 
    plt.plot(time,series)
    plt.show()

    my_series = TimeSeries(series=series, step=dt)
    
    chrystals_dft_freqs, chrystals_dft_amplitudes = my_series.chrystals_dft()
    
    N_power2 = int(2**(np.floor(np.log2(len(series)))))
    chrystals_fft_amplitudes = 2.0/N_power2*np.abs(chrystals_fft(series=series[:N_power2]))[1:N_power2//2]
    chrystals_fft_freqs = np.arange(1,N_power2//2)/N_power2/dt
    
    from scipy.fft import fft, fftfreq
    frequencies = fftfreq(nt,dt)[1:nt//2]
    amplitudes = 2.0/nt*np.abs(fft(series)[1:nt//2])

    from matplotlib import pyplot as plt
    plt.plot(frequencies, amplitudes, label="scipy")
    plt.plot(chrystals_fft_freqs, chrystals_fft_amplitudes, 'o', label="fft")
    plt.plot(chrystals_dft_freqs, chrystals_dft_amplitudes, 'x', label="dft")
    plt.legend()
    plt.show()
