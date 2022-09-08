import numpy as np
from sklearn.neighbors import KernelDensity

def estimate_univariate_kde(posteriorsamples, name):
    KDE = []
    for i in np.arange(0, len(posteriorsamples[1, :])):
        print(i)
        log_estimated_parameters = posteriorsamples[:, i]
        silverman_bandwidth = (4/3)**(1/5)*np.std(log_estimated_parameters)*len(log_estimated_parameters)**(-1/5)
        kdei = KernelDensity(bandwidth=silverman_bandwidth, kernel='gaussian').fit(log_estimated_parameters.reshape(-1, 1))
        KDE.append(kdei)
    np.save('kdes_{}.npy'.format(name), KDE)