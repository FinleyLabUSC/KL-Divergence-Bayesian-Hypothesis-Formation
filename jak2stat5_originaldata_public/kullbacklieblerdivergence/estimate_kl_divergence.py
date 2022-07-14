import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.stats import anderson
import statsmodels.api as sm
import pandas as pd
import seaborn as sns
from sklearn.cluster import KMeans

def estimate_kl_divergence_posterior_prior(posteriorsamples_all, prior_std_dev, KDE):
    prior_mu = np.log(np.load('kullbacklieblerdivergence/prior_estimatedvariables.npy'))
    variable_dictionary = np.load('kullbacklieblerdivergence/dictionary_estimatedvariables.npy')
    posteriorsamples = posteriorsamples_all

    KL = []
    for i in np.arange(0, len(KDE)):
        print(i)
        kde = KDE[i]
        pri_mu = prior_mu[i]
        samples = np.log(posteriorsamples[:, i])
        div = []
        for s in samples:
            log_prior = norm.logpdf(s, loc=pri_mu, scale=prior_std_dev)
            # note - we use the normal distribution for evaluating prior probability,
            # because taking a lognormal distribution into logscale makes it the normal distribution.
            # normal distribution is better than long-tailed lognormal for kde methods
            log_posterior = kde.score_samples(s.reshape(1, -1))
            div.append(log_posterior - log_prior)
        KLi = np.mean(np.array(div))
        KL.append(KLi)
    np.save('kullbacklieblerdivergence/kullback_leibler_divergence_posterior_prior.npy', KL)

def KL_div_scatter_plot_clustered(kl_div, fontsize):
    variable_dictionary = np.load('kullbacklieblerdivergence/dictionary_estimatedvariables.npy')
    kmeans = KMeans(2).fit(kl_div.reshape(-1, 1))
    labels = kmeans.labels_
    data = np.concatenate([kl_div.reshape(33, 1), variable_dictionary.reshape(33, 1), labels.reshape(33, 1)], axis=1)
    DF = pd.DataFrame(data=data, columns=['KL', 'Names', 'Cluster'])
    DF_sort = DF.sort_values(by='KL')
    sorted_labels = np.array(list(map(float, DF_sort['Cluster'])))
    cluster_one = np.argwhere(sorted_labels==1)
    cluster_two = np.argwhere(sorted_labels==0)
    all = np.array(list(map(float, DF_sort['KL'])))
    dictionary = np.array(list(DF_sort['Names']))
    plt.figure(figsize=(10, 7))
    plt.plot(np.arange(1, 34), np.around(all, decimals=3), 'o', color='black', markersize=10)
    plt.plot(np.arange(1, len(cluster_two)+1), np.around(all[cluster_two], decimals=3), 'o', color='slategrey', markersize=10)
    plt.xticks(ticks=np.arange(1, 34), labels=dictionary, rotation='vertical')
    plt.ylabel('Kullback-Leibler Divergence', fontsize=fontsize)
    plt.xlabel('Parameters', fontsize=fontsize)
    plt.tight_layout()
    plt.savefig('kullbacklieblerdivergence/kullback_leibler_divergence_posterior_prior_cluster.png', dpi=300)
    plt.close()

