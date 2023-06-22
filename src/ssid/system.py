from ssid import markov, realize, modal

def system(input, output, method="srim", **options):

    if method == "okid-era":
        Y = markov.okid(input, output, **options)
        realization = realize.era(Y, **options)

    if method == "okid-era-dc":
        Y = markov.okid(input, output, **options)
        realization = realize.era_dc(Y, **options)
    
    if method == "srim":
        realization = realize.srim(input, output, **options)

    return realization

# def spectrum():
    # return

def modes(input, output, dt, method="srim", **options):
    
    realization = system(input, output, method, **options)

    return modal.system_modes(realization, dt)

