import numpy as np
from scipy.stats import t
from scipy.stats import norm
import math

def log_likelihood_baseline_normal_inverse_gamma_variance(x):
    bcl = x[1:8]
    pjak = x[9:15]
    pstata = x[17:23]
    pstatb = x[25:31]
    transloca = x[32:37]
    translocb = x[37:46]
    exp_data = np.load('posteriorapproximation/data_experimental.npy').flatten()
    exp_mu_bcl = exp_data[26:33]
    exp_mu_pjak = exp_data[33:39]
    exp_mu_pstata = exp_data[0:6]
    exp_mu_pstatb = exp_data[6:12]
    exp_mu_transloca = exp_data[12:17]
    exp_mu_translocb = exp_data[17:26]
    predictions = np.concatenate((bcl, pjak, pstata, pstatb, transloca, translocb))
    data = np.concatenate((exp_mu_bcl, exp_mu_pjak, exp_mu_pstata, exp_mu_pstatb, exp_mu_transloca, exp_mu_translocb))
    error = np.sum((np.divide(np.subtract(predictions, data), data))**2)
    alpha = 2
    beta = 1e-3
    datasize = len(data)
    ll = - (datasize / 2 + alpha) * math.log(error / 2 + beta)
    print('error', error)
    return ll

def log_likelihood_baseline_and_timepoint(x):
    bcl = x[1:8]
    pjak = x[9:16]
    pstata = x[17:24]
    pstatb = x[25:32]
    transloca = x[32:37]
    translocb = x[37:46]
    exp_data = np.load('posteriorapproximation/data_experimental.npy').flatten()
    exp_mu_pjak_12 = np.load('1. Data/exp_pjak_mu.npy')[-1]
    exp_mu_pstata_12 = np.load('1. Data/exp_pstata_mu.npy')[-1]
    exp_mu_pstatb_12 = np.load('1. Data/exp_pstatb_mu.npy')[-1]
    exp_mu_bcl = exp_data[26:33]
    exp_mu_pjak = np.append(exp_data[33:39], exp_mu_pjak_12)
    exp_mu_pstata = np.append(exp_data[0:6], exp_mu_pstata_12)
    exp_mu_pstatb = np.append(exp_data[6:12], exp_mu_pstatb_12)
    exp_mu_transloca = exp_data[12:17]
    exp_mu_translocb = exp_data[17:26]
    predictions = np.concatenate((bcl, pjak, pstata, pstatb, transloca, translocb))
    data = np.concatenate((exp_mu_bcl, exp_mu_pjak, exp_mu_pstata, exp_mu_pstatb, exp_mu_transloca, exp_mu_translocb))
    error = np.sum((np.divide(np.subtract(predictions, data), data))**2)
    alpha = 2
    beta = 1e-3
    datasize = len(data)
    ll = - (datasize / 2 + alpha) * math.log(error / 2 + beta)
    print('error', error)
    return ll



