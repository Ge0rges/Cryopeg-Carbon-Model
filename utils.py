def convert_micromolar_carbon_to_fgC_per_ml (micromolar):
    return (micromolar * 12.011 * 10**-6 * 10**15)/1000

def convert_fgC_per_ml_to_micromolar_carbon(fgC_per_ml):
    return (fgC_per_ml * 1000)/(12.011 * 10**-6 * 10**15)
