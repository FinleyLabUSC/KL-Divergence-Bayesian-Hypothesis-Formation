import numpy as np
from scipy.signal import find_peaks
from sklearn.neighbors import KernelDensity
import matplotlib.pyplot as plt
import seaborn as sns

#code modes in marginal posteriors based on their univariate KDEs

#KDEs
posteriorsamples = np.load('generateinputs/all_chains_baseline_proposal 0.1_thinnedestimatedvariables.npy')
KDE = []
for i in np.arange(0, len(posteriorsamples[1, :])):
    print(i)
    log_estimated_parameters = np.log(posteriorsamples[:, i])
    silverman_bandwidth = (4 / 3) ** (1 / 5) * np.std(log_estimated_parameters) * len(log_estimated_parameters) ** (
                -1 / 5)
    kdei = KernelDensity(bandwidth=silverman_bandwidth, kernel='gaussian').fit(log_estimated_parameters.reshape(-1, 1))
    KDE.append(kdei)

#plot bimodalities
variable_dictionary = np.load('kullbacklieblerdivergence/dictionary_estimatedvariables.npy')
fig = plt.figure(figsize=(8, 11))
for i in np.arange(0, 33):
    print(i)
    ax = fig.add_subplot(7, 5, i + 1)
    plt.subplots_adjust(bottom=0.05, top=0.95, left=0.05, right=0.99)
    estimated_parameters = np.log(posteriorsamples[:, i])
    sorted_estimated_parameters = np.sort(estimated_parameters)
    KDEi = KDE[i]
    scores = np.exp(KDEi.score_samples(sorted_estimated_parameters.reshape(-1, 1)))
    peaks, properties = find_peaks(scores)
    peak_y = scores[peaks]
    peak_x = sorted_estimated_parameters[peaks]
    mode_count = len(peaks)
    # sns.histplot(data=estimated_parameters, color='slategrey', ax=ax)
    # sns.histplot(data=scores, color='dodgerblue', alpha=0.3, ax=ax)
    plt.plot(sorted_estimated_parameters, scores, color='dodgerblue')
    plt.plot(peak_x, peak_y, '*', color='firebrick')
    plt.xlabel('ln({})'.format(variable_dictionary[i]))
    plt.ylabel('')
    plt.yticks([0], '')
    # plt.ylim([0, 12000])
plt.tight_layout()
plt.savefig('multimodality_test/all_modes.png', dpi=300)
plt.close(fig)