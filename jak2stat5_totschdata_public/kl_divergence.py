import numpy as np
from scipy.stats import norm
from scipy.stats import uniform


def estimate_kl_divergence_posterior_loguninform_prior_totsch(posteriorsamples_all, prior_lower, prior_upper, KDE, name):
    posteriorsamples = posteriorsamples_all[0:-1:int(np.rint(len(posteriorsamples_all)/10000))]
    div = []
    for s in posteriorsamples:
        log_prior = uniform.logpdf(s, loc=prior_lower, scale=prior_upper - prior_lower)
        log_posterior = KDE.score_samples(s.reshape(1, -1))
        div.append(log_posterior - log_prior)
    KLi = np.mean(np.array(div))
    np.save('kullback_leibler_divergence_{}.npy'.format(name), KLi)
