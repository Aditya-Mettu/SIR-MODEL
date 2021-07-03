import numpy as np
from scipy.optimize import curve_fit

#
# 'exp1' for exponential one model.
# 'logistic' for logistic model.
# 'logmdl' for 'for coefficients' model.
#

def funcExp1(x, a, b):
    return a*np.exp(b*x)

def funcLogistic(x, a1, a2, a3):
    return a1/(1 + a2*np.exp(-1*a3*x))

def funcLogmdl(x, a0, a1, a2):
    return a0/(1+a1*np.exp(-1*a2*x))

def fit(texp, Cexp, name):
    
    popt, pcov = None, None
    if name == 'exp1':
        popt, pcov = curve_fit(funcExp1, texp, Cexp)
    elif name == 'logistics':
        popt, pcov = curve_fit(funcLogistic, texp, Cexp)
    elif name == 'logmdl':
        popt, pcov = curve_fit(funcLogmdl, texp, Cexp)

    return popt