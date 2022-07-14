import numpy as np
import math

def generate_lognormal(parameters, std):
    samples = []
    for item in parameters:
        if item <= 1e-10:
            item = 1e-10
        item = math.log(item)
        sample = np.random.lognormal(item, std)
        samples.append(sample)
    return samples

def evaluate_lognormal_function(parameters, mu, std):
    mu = np.log(mu)
    p = -(1/(2*std**2))*np.sum((np.log(parameters)-mu)**2)
    return p

