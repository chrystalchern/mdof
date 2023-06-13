from ssid import markov, realize, modes

def system(input, output, method="srim", **options):

    if method == "okid-era":
        Y = markov.okid(input, output, **options)
        realilzation = realize.era(Y, **options)

    if method == "okid-era-dc":
        Y = markov.okid(input, output, **options)
        realization = realize.era_dc(Y, **options)
    
    if method == "srim":
        realization = realize.srim(input, output, **options)

    return realization

def spectrum(input, output, dt, method="srim", **options):
    
    realization = system(input, output, method, **options)

    return modes.modes(realization, dt)