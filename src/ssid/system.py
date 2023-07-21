from ssid import markov, realize, modal
from .numerics import decimate

def system(inputs, outputs, method="srim", decimation=None, **options):

    if method not in {
        "srim",
        "okid-era",
        "okid-era-dc"
    }: raise ValueError(f"Unknown method {method}")

    if decimation is not None:
        inputs = decimate(inputs, decimation=decimation)
        outputs = decimate(outputs, decimation=decimation)

    if method == "okid-era":
        Y = markov.okid(inputs, outputs, **options)
        realization = realize.era(Y, **options)

    if method == "okid-era-dc":
        Y = markov.okid(inputs, outputs, **options)
        realization = realize.era_dc(Y, **options)
    
    if method == "srim":
        realization = realize.srim(inputs, outputs, **options)

    return realization

def modes(inputs, outputs, dt, method="srim", decimation=None, **options):

    realization = system(inputs, outputs, method, decimation, **options)

    return modal.system_modes(realization, dt, decimation)

