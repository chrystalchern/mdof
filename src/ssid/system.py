from ssid import markov, realize

def system(method="srim", *args, **options):

    if method == "okid-era":
        Y = markov.okid(*args, **options)
        A,B,C,D = realize.era(Y, **options)

    if method == "okid-era-dc":
        Y = markov.okid(*args, **options)
        A,B,C,D = realize.era_dc(Y, **options)
    
    if method == "srim":
        A,B,C,D = realize.srim(*args, **options)

    return A,B,C,D