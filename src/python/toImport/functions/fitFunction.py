import numpy as np
from scipy.optimize import curve_fit

def funcExp1(x, a, b):
    return a*np.exp(b*x)

def fit(texp, Cexp, name):
    
    popt, pcov = None, None
    if name == 'exp1':
        popt, pcov = curve_fit(funcExp1, texp, Cexp)

    return popt