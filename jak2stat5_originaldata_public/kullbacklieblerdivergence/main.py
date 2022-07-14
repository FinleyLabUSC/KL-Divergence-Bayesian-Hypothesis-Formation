from kernel_density_estimate import estimate_univariate_kde
from estimate_kl_divergence import estimate_kl_divergence_posterior_prior
from estimate_kl_divergence import KL_div_scatter_plot_clustered
import numpy as np
from sklearn.cluster import KMeans
from sklearn import metrics
from sklearn.metrics import pairwise_distances
import matplotlib.pyplot as plt

#estimate KDEs - also generates figure checking KDE samples against original samples,
# displayed in Supplementary Materials section 2
posterior_samples = np.load('generateinputs/all_chains_baseline_proposal 0.1_thinnedestimatedvariables.npy')
estimate_univariate_kde(posterior_samples)

#estimate KL Divergence
priorstddev = np.sqrt(2)
KDE = np.load('kullbacklieblerdivergence/posterior_KDEs.npy', allow_pickle=True)
estimate_kl_divergence_posterior_prior(posterior_samples, priorstddev, KDE)

#cluster KL divergence
KL_D = np.load('kullbacklieblerdivergence/kullback_leibler_divergence_posterior_prior.npy')
fontsize = 18

#test cluster - print out names of two groups
kmeans2 = KMeans(2).fit(KL_D.reshape(-1, 1))
labels = kmeans2.labels_
dictionary = np.load('kullbacklieblerdivergence/dictionary_estimatedvariables.npy')
print("group 1: \n",  dictionary[np.argwhere(labels==1)])
print("group 0: \n",  dictionary[np.argwhere(labels==0)])
#note - group number assignment is arbitrary

#cluster and generate figure 3A
KL_div_scatter_plot_clustered(KL_D, fontsize)

#generate elbow plot, silhouette score to determine number of clusters (Supplementary Figure S2)
labels = []
WCSS = []
SS = []

for i in range(2, 8):
    #note, must run for loop seperately for WCSS when k=1, because silhouette score does not apply for k=1
    kmeans = KMeans(i).fit(KL_D.reshape(-1, 1))
    WCSS.append(kmeans.inertia_)
    labels.append(kmeans.labels_)
    silhouette_score = metrics.silhouette_score(KL_D.reshape(-1, 1), kmeans.labels_, metric='euclidean')
    SS.append(silhouette_score)

WCSS = np.load('kullbacklieblerdivergence/WCSS.npy')
plt.plot(np.arange(1, 8), WCSS, 'o-', color='slategrey')
plt.plot(2, WCSS[1],'o', color='firebrick')
plt.xlabel('Number of Clusters')
plt.ylabel('Within Cluster Sum of Squares Criterion')
for i in np.arange(0, 6):
    plt.annotate(np.around(SS[i], 2), xy=(i + 2 + 0.1, WCSS[i+1] + 0.1), color='dodgerblue')
plt.xlim([0, 8])
plt.savefig('kullbacklieblerdivergence/elbow_plot.png', dpi=300)
plt.show()


