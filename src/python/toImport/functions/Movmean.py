import numpy as np 

def movmean(a, n=3) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return np.concatenate((a[:n-1], ret[n - 1:] / n), axis=0)