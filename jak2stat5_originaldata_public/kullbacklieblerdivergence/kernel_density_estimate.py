import numpy as np
from sklearn.neighbors import KernelDensity
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.signal import find_peaks

def estimate_univariate_kde(posteriorsamples):
    KDE = []
    for i in np.arange(0, len(posteriorsamples[1, :])):
        print(i)
        log_estimated_parameters = np.log(posteriorsamples[:, i])
        silverman_bandwidth = (4/3)**(1/5)*np.std(log_estimated_parameters)*len(log_estimated_parameters)**(-1/5)
        kdei = KernelDensity(bandwidth=silverman_bandwidth, kernel='gaussian').fit(log_estimated_parameters.reshape(-1, 1))
        KDE.append(kdei)
    np.save('kullbacklieblerdivergence/posterior_KDEs.npy', KDE)

    variable_dictionary = np.load('kullbacklieblerdivergence/dictionary_estimatedvariables.npy')

    fig = plt.figure(figsize=(13, 10))
    for i in np.arange(0, 11):
        print(i)
        ax = fig.add_subplot(3, 4, i + 1)
        plt.subplots_adjust(bottom=0.05, top=0.95, left=0.05, right=0.99)
        estimated_parameters = np.log(posteriorsamples[:, i])
        KDEi = KDE[i]
        samples = KDEi.sample(len(posteriorsamples[:, i]))
        sns.histplot(data=estimated_parameters, color='slategrey', ax=ax)
        sns.histplot(data=samples, color='dodgerblue', alpha=0.3, ax=ax)
        if i==0:
            plt.legend(['true posterior samples', 'kernel posterior samples'])
        else:
            plt.legend([])
        #if i != 0 and i != 4 and i != 8:
            #plt.yticks([0, 12000], '')
        plt.xlabel('ln({})'.format(variable_dictionary[i]))
        plt.ylabel('')
        #plt.ylim([0, 12000])
    plt.savefig('kullbacklieblerdivergence/check_kde/group1.png', dpi=300)
    plt.close(fig)

    fig = plt.figure(figsize=(13, 10))
    for i in np.arange(11, 22):
        print(i)
        ax = fig.add_subplot(3, 4, i-11 + 1)
        plt.subplots_adjust(bottom=0.05, top=0.95, left=0.05, right=0.99)
        estimated_parameters = np.log(posteriorsamples[:, i])
        KDEi = KDE[i]
        samples = KDEi.sample(len(posteriorsamples[:, i]))
        sns.histplot(data=estimated_parameters, color='slategrey', ax=ax)
        sns.histplot(data=samples, color='dodgerblue', alpha=0.3, ax=ax)
        if i==11:
            plt.legend(['true posterior samples', 'kernel posterior samples'])
        else:
            plt.legend([])
        #if i != 0+11 and i != 4+11 and i != 8+11:
            #plt.yticks([0, 12000], '')
        plt.xlabel('ln({})'.format(variable_dictionary[i]))
        plt.ylabel('')
        #plt.ylim([0, 12000])
    plt.savefig('kullbacklieblerdivergence/check_kde/group2.png', dpi=300)
    plt.close(fig)

    fig = plt.figure(figsize=(13, 10))
    for i in np.arange(22, 33):
        print(i)
        ax = fig.add_subplot(3, 4, i-22 + 1)
        plt.subplots_adjust(bottom=0.05, top=0.95, left=0.05, right=0.99)
        estimated_parameters = np.log(posteriorsamples[:, i])
        KDEi = KDE[i]
        samples = KDEi.sample(len(posteriorsamples[:, i]))
        sns.histplot(data=estimated_parameters, color='slategrey', ax=ax)
        sns.histplot(data=samples, color='dodgerblue', alpha=0.3, ax=ax)
        if i==22:
            plt.legend(['true posterior samples', 'kernel posterior samples'])
        else:
            plt.legend([])
        #if i != 0+22 and i != 4+22 and i != 8+22:
            #plt.yticks([0, 13000], '')
        plt.xlabel('ln({})'.format(variable_dictionary[i]))
        plt.ylabel('')
        #plt.ylim([0, 12000])
    plt.savefig('kullbacklieblerdivergence/check_kde/group3.png', dpi=300)
    plt.close(fig)

    fig = plt.figure(figsize=(13, 10))
    for i in np.arange(0, 11):
        print(i)
        ax = fig.add_subplot(3, 4, i + 1)
        plt.subplots_adjust(bottom=0.05, top=0.95, left=0.05, right=0.99)
        estimated_parameters = np.log(posteriorsamples[:, i])
        sorted_estimated_parameters = np.sort(estimated_parameters)
        KDEi = KDE[i]
        scores = np.exp(KDEi.score_samples(sorted_estimated_parameters.reshape(-1, 1)))
        peaks, properties = find_peaks(scores)
        peak_y = scores[peaks]
        peak_x = sorted_estimated_parameters[peaks]
        mode_count = len(peaks)
        #sns.histplot(data=estimated_parameters, color='slategrey', ax=ax)
        #sns.histplot(data=scores, color='dodgerblue', alpha=0.3, ax=ax)
        plt.plot(sorted_estimated_parameters, scores, color='dodgerblue')
        plt.plot(peak_x,peak_y,'*', color='firebrick')
        if i==0:
            plt.legend(['kernel density', 'mode count: {}'.format(mode_count)], loc='lower center')
        else:
            plt.legend([])
        plt.xlabel('ln({})'.format(variable_dictionary[i]))
        plt.ylabel('')
        #plt.ylim([0, 12000])
    plt.savefig('kullbacklieblerdivergence/mode_count/group1_modes.png', dpi=300)
    plt.close(fig)

    fig = plt.figure(figsize=(13, 10))
    for i in np.arange(11, 22):
        print(i)
        ax = fig.add_subplot(3, 4, i-11 + 1)
        plt.subplots_adjust(bottom=0.05, top=0.95, left=0.05, right=0.99)
        estimated_parameters = np.log(posteriorsamples[:, i])
        sorted_estimated_parameters = np.sort(estimated_parameters)
        KDEi = KDE[i]
        scores = np.exp(KDEi.score_samples(sorted_estimated_parameters.reshape(-1, 1)))
        peaks, properties = find_peaks(scores)
        peak_y = scores[peaks]
        peak_x = sorted_estimated_parameters[peaks]
        mode_count = len(peaks)
        #sns.histplot(data=estimated_parameters, color='slategrey', ax=ax)
        #sns.histplot(data=samples, color='dodgerblue', alpha=0.3, ax=ax)
        plt.plot(sorted_estimated_parameters, scores, color='dodgerblue')
        plt.plot(peak_x,peak_y,'*', color='firebrick')
        if i==11:
            plt.legend(['kernel density', 'mode count: {}'.format(mode_count)], loc='lower center')
        else:
            plt.legend([])
        plt.xlabel('ln({})'.format(variable_dictionary[i]))
        plt.ylabel('')
        #plt.ylim([0, 12000])
    plt.savefig('kullbacklieblerdivergence/mode_count/group2_modes.png', dpi=300)
    plt.close(fig)

    fig = plt.figure(figsize=(13, 10))
    for i in np.arange(22, 33):
        print(i)
        ax = fig.add_subplot(3, 4, i-22 + 1)
        plt.subplots_adjust(bottom=0.05, top=0.95, left=0.05, right=0.99)
        estimated_parameters = np.log(posteriorsamples[:, i])
        sorted_estimated_parameters = np.sort(estimated_parameters)
        KDEi = KDE[i]
        scores = np.exp(KDEi.score_samples(sorted_estimated_parameters.reshape(-1, 1)))
        peaks, properties = find_peaks(scores)
        peak_y = scores[peaks]
        peak_x = sorted_estimated_parameters[peaks]
        mode_count = len(peaks)
        #sns.histplot(data=estimated_parameters, color='slategrey', ax=ax)
        #sns.histplot(data=samples, color='dodgerblue', alpha=0.3, ax=ax)
        plt.plot(sorted_estimated_parameters, scores, color='dodgerblue')
        plt.plot(peak_x,peak_y,'*', color='firebrick')
        if i==22:
            plt.legend(['kernel density', 'mode count: {}'.format(mode_count)], loc='lower center')
        else:
            plt.legend([])
        plt.xlabel('ln({})'.format(variable_dictionary[i]))
        plt.ylabel('')
        #plt.ylim([0, 12000])
    plt.savefig('kullbacklieblerdivergence/mode_count/group3_modes.png', dpi=300)
    plt.close(fig)

    fig = plt.figure(figsize=(13, 10))
    for i in np.arange(0, 11):
        print(i)
        ax = fig.add_subplot(3, 4, i + 1)
        plt.subplots_adjust(bottom=0.05, top=0.95, left=0.05, right=0.99)
        estimated_parameters = np.log(posteriorsamples[:, i])
        sorted_estimated_parameters = np.sort(estimated_parameters)
        KDEi = KDE[i]
        scores = np.exp(KDEi.score_samples(sorted_estimated_parameters.reshape(-1, 1)))
        #sns.histplot(data=estimated_parameters, color='slategrey', ax=ax)
        #sns.histplot(data=scores, color='dodgerblue', alpha=0.3, ax=ax)
        plt.plot(sorted_estimated_parameters, scores, color='dodgerblue')
        if i==0:
            plt.legend(['kernel density'], loc='lower center')
        else:
            plt.legend([])
        plt.xlabel('ln({})'.format(variable_dictionary[i]))
        plt.ylabel('')
        #plt.ylim([0, 12000])
    plt.savefig('kullbacklieblerdivergence/mode_count/group1_kde.png', dpi=300)
    plt.close(fig)

    fig = plt.figure(figsize=(13, 10))
    for i in np.arange(11, 22):
        print(i)
        ax = fig.add_subplot(3, 4, i-11 + 1)
        plt.subplots_adjust(bottom=0.05, top=0.95, left=0.05, right=0.99)
        estimated_parameters = np.log(posteriorsamples[:, i])
        sorted_estimated_parameters = np.sort(estimated_parameters)
        KDEi = KDE[i]
        scores = np.exp(KDEi.score_samples(sorted_estimated_parameters.reshape(-1, 1)))
        #sns.histplot(data=estimated_parameters, color='slategrey', ax=ax)
        #sns.histplot(data=samples, color='dodgerblue', alpha=0.3, ax=ax)
        plt.plot(sorted_estimated_parameters, scores, color='dodgerblue')
        if i==11:
            plt.legend(['kernel density'], loc='lower center')
        else:
            plt.legend([])
        plt.xlabel('ln({})'.format(variable_dictionary[i]))
        plt.ylabel('')
        #plt.ylim([0, 12000])
    plt.savefig('kullbacklieblerdivergence/mode_count/group2_kde.png', dpi=300)
    plt.close(fig)

    fig = plt.figure(figsize=(13, 10))
    for i in np.arange(22, 33):
        print(i)
        ax = fig.add_subplot(3, 4, i-22 + 1)
        plt.subplots_adjust(bottom=0.05, top=0.95, left=0.05, right=0.99)
        estimated_parameters = np.log(posteriorsamples[:, i])
        sorted_estimated_parameters = np.sort(estimated_parameters)
        KDEi = KDE[i]
        scores = np.exp(KDEi.score_samples(sorted_estimated_parameters.reshape(-1, 1)))
        #sns.histplot(data=estimated_parameters, color='slategrey', ax=ax)
        #sns.histplot(data=samples, color='dodgerblue', alpha=0.3, ax=ax)
        plt.plot(sorted_estimated_parameters, scores, color='dodgerblue')
        if i==22:
            plt.legend(['kernel density'], loc='lower center')
        else:
            plt.legend([])
        plt.xlabel('ln({})'.format(variable_dictionary[i]))
        plt.ylabel('')
        #plt.ylim([0, 12000])
    plt.savefig('kullbacklieblerdivergence/mode_count/group3_kde.png', dpi=300)
    plt.close(fig)

