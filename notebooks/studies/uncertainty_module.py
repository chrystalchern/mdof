import numpy as np
from datascience import *

def quantify_noise(series, expected_series):
    # Temporary stand-in for quantifying the noise.  We'll need a better function later.
    return np.mean(series - expected_series)/np.max(expected_series)


def determine_uncertainty(inputs: np.array, outputs: np.array, dt):
    # Calculations for expected input and output
    p, nt = outputs.shape
    q, nt = inputs.shape
    time = np.arange(0,nt*dt,dt)
    expected_output = -0.172*np.exp(-0.05*time)*np.cos(3.16*time) + 0.025*np.exp(-0.05*time)*np.sin(3.16*time) + 0.172*np.exp(-time)*np.cos(np.pi*time) + 0.027*np.exp(-time)*np.sin(np.pi*time)
    expected_input = np.exp(-time)*np.sin(np.pi*time)

    # Calculations for amount of noise, comparing measured vibration to expected vibration
    noise_output =  quantify_noise(outputs, expected_output)         
    noise_input =   quantify_noise(inputs, expected_input)

    # Give a level of uncertainty.

    # Method 1: Give hard-coded numbers.
    if noise_output >= 1 and noise_input >= 1:   # A lot of output noise AND a lot of input noise
        level_of_certainty = ...
    elif noise_output >= 1 and noise_input < 1:   # A lot of output noise but not a lot of input noise
        ...
    elif ...:
        ...
    else:
        ...

    # Method 2: Calculate a level of certainty that decreases as noise increases
    level_of_certainty = -noise_output + ...

    return level_of_certainty


if __name__ == '__main__':
    time = np.array([np.arange(0,50,0.1)])
    inputs = np.exp(-time)*np.sin(np.pi*time) + np.random.normal(0,0.01,500)
    outputs = -0.172*np.exp(-0.05*time)*np.cos(3.16*time) + 0.025*np.exp(-0.05*time)*np.sin(3.16*time) + 0.172*np.exp(-time)*np.cos(np.pi*time) + 0.027*np.exp(-time)*np.sin(np.pi*time) + np.random.normal(0,0.0001,500)
    determine_uncertainty(inputs, outputs, 0.1)
