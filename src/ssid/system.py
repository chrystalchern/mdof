from ssid import markov, realize, modal

def system(inputs, outputs, method="srim", **options):

    if method == "okid-era":
        Y = markov.okid(inputs, outputs, **options)
        realization = realize.era(Y, **options)

    if method == "okid-era-dc":
        Y = markov.okid(inputs, outputs, **options)
        realization = realize.era_dc(Y, **options)
    
    if method == "srim":
        realization = realize.srim(inputs, outputs, **options)

    return realization

# def spectrum():
    # return

def modes(inputs, outputs, dt, method="srim", **options):
    
    realization = system(inputs, outputs, method, **options)

    return modal.system_modes(realization, dt)

