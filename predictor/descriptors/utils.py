import numpy as np

def impute(X, medians):
    for i in range(X.shape[1]):
        vals = X[:,i]
        mask = np.isnan(vals)
        vals[mask] = medians[i]
        X[:,i] = vals
    return X

def normalize(X, mus, sigmas):
    for i in range(X.shape[1]):
        vals = X[:,i]
        sigma = sigmas[i]
        if sigma == 0:
            vals = np.zeros(len(vals))
        else:
            vals = (vals - mus[i])/sigmas[i]
        X[:,i] = vals
    return X
