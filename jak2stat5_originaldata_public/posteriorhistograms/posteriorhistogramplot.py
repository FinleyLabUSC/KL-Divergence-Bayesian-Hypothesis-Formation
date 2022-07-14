import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def histogramsposterior_logscale(posteriorsamples, path):
    # histogram with 95th percentile of prior
    variable_dictionary = np.load('posteriorhistograms/dictionary_estimatedvariables.npy')
    k24_index = np.where(variable_dictionary == 'k24')[0][0]

    plt.figure(figsize=(13, 10))
    for i in np.arange(0, 11):
        parameter = np.log(posteriorsamples[:, i])
        plt.subplot(3, 4, i + 1)
        plt.subplots_adjust(bottom=0.05, top=0.95, left=0.05, right=0.99)
        sns.histplot(data=parameter, color='slategrey', kde=True)
        plt.xlabel('ln({})'.format(variable_dictionary[i]), fontsize=14)
        if i != 0 and i != 4 and i != 8:
            plt.yticks([0, 25000], '')
        plt.ylabel('', fontsize=14)
        plt.ylim([0, 25000])
    plt.savefig('{}/histogram_group1_posterior.png'.format(path), dpi=150)

    plt.figure(figsize=(13, 10))
    for i in np.arange(11, 22):
        parameter = np.log(posteriorsamples[:, i])
        plt.subplot(3, 4, i - 11 + 1)
        plt.subplots_adjust(bottom=0.05, top=0.95, left=0.05, right=0.99)
        if i == k24_index:
            sns.histplot(data=parameter, color='black', kde=True)
        else:
            sns.histplot(data=parameter, color='slategrey', kde=True)
        plt.xlabel('ln({})'.format(variable_dictionary[i]), fontsize=14)
        if i != 0+11 and i != 4+11 and i != 8+11:
            plt.yticks([0, 25000], '')
        plt.ylabel('', fontsize=14)
        plt.ylim([0, 25000])
    plt.savefig('{}/histogram_group2_posterior.png'.format(path), dpi=150)

    plt.figure(figsize=(13, 10))
    for i in np.arange(22, 33):
        parameter = np.log(posteriorsamples[:, i])
        plt.subplot(3, 4, i - 22 + 1)
        plt.subplots_adjust(bottom=0.05, top=0.95, left=0.05, right=0.99)
        sns.histplot(data=parameter, color='slategrey', kde=True)
        plt.xlabel('ln({})'.format(variable_dictionary[i]), fontsize=14)
        if i != 0+22 and i != 4+22 and i != 8+22:
            plt.yticks([0, 25000], '')
        plt.ylabel('', fontsize=14)
        plt.ylim([0, 25000])
    plt.savefig('{}/histogram_group3_posterior.png'.format(path), dpi=150)

def histogramspriorposterior_logscale(posteriorsamples, prior_std_dev, path):
    # histogram with 95th percentile of prior
    prior_estimatedvariables = np.load('posteriorhistograms/prior_estimatedvariables.npy')
    variable_dictionary = np.load('posteriorhistograms/dictionary_estimatedvariables.npy')
    std_prior = prior_std_dev
    k24_index = np.where(variable_dictionary == 'k24')[0][0]

    plt.figure(figsize=(13, 10))
    for i in np.arange(0, 11):
        parameter = np.log(posteriorsamples[:, i])
        mu = np.log(prior_estimatedvariables[i])
        prior = np.log(np.random.lognormal(mean=mu, sigma=std_prior, size=len(parameter)))
        upper_percentile = mu + 2 * std_prior
        lower_percentile = mu - 2 * std_prior
        plt.subplot(3, 4, i + 1)
        plt.subplots_adjust(bottom=0.05, top=0.95, left=0.05, right=0.99)
        sns.histplot(data=parameter, color='slategrey')
        sns.histplot(data=prior, color='slategrey', kde=True)
        plt.vlines(x=upper_percentile, ymin=0, ymax=25000, color='darkturquoise', linestyles='dashed')
        plt.vlines(x=lower_percentile, ymin=0, ymax=25000, color='darkturquoise', linestyles='dashed')
        plt.xlabel('ln({})'.format(variable_dictionary[i]), fontsize=14)
        if i != 0 and i != 4 and i != 8:
            plt.yticks([0, 25000], '')
        plt.ylabel('', fontsize=14)
        plt.ylim([0, 25000])
    plt.savefig('{}/histogram_group1_prior&posterior.png'.format(path))

    plt.figure(figsize=(13, 10))
    for i in np.arange(11, 22):
        parameter = np.log(posteriorsamples[:, i])
        mu = np.log(prior_estimatedvariables[i])
        prior = np.log(np.random.lognormal(mean=mu, sigma=std_prior, size=len(parameter)))
        upper_percentile = mu + 2 * std_prior
        lower_percentile = mu - 2 * std_prior
        plt.subplot(3, 4, i - 11 + 1)
        plt.subplots_adjust(bottom=0.05, top=0.95, left=0.05, right=0.99)
        if i == k24_index:
            sns.histplot(data=parameter, color='black')
            sns.histplot(data=prior, color='black', kde=True)
        else:
            sns.histplot(data=parameter, color='slategrey')
            sns.histplot(data=prior, color='slategrey', kde=True)
        plt.vlines(x=upper_percentile, ymin=0, ymax=25000, color='darkturquoise', linestyles='dashed')
        plt.vlines(x=lower_percentile, ymin=0, ymax=25000, color='darkturquoise', linestyles='dashed')
        plt.xlabel('ln({})'.format(variable_dictionary[i]), fontsize=14)
        if i != 0 and i != 4 and i != 8:
            plt.yticks([0, 25000], '')
        plt.ylabel('', fontsize=14)
        plt.ylim([0, 25000])
    plt.savefig('{}/histogram_group2_prior&posterior.png'.format(path))

    plt.figure(figsize=(13, 10))
    for i in np.arange(22, 33):
        parameter = np.log(posteriorsamples[:, i])
        mu = np.log(prior_estimatedvariables[i])
        prior = np.log(np.random.lognormal(mean=mu, sigma=std_prior, size=len(parameter)))
        upper_percentile = mu + 2 * std_prior
        lower_percentile = mu - 2 * std_prior
        plt.subplot(3, 4, i - 22 + 1)
        plt.subplots_adjust(bottom=0.05, top=0.95, left=0.05, right=0.99)
        sns.histplot(data=parameter, color='slategrey')
        sns.histplot(data=prior, color='slategrey', kde=True)
        plt.vlines(x=upper_percentile, ymin=0, ymax=25000, color='darkturquoise', linestyles='dashed')
        plt.vlines(x=lower_percentile, ymin=0, ymax=25000, color='darkturquoise', linestyles='dashed')
        plt.xlabel('ln({})'.format(variable_dictionary[i]), fontsize=14)
        if i != 0 and i != 4 and i != 8:
            plt.yticks([0, 25000], '')
        plt.ylabel('', fontsize=14)
        plt.ylim([0, 25000])
    plt.savefig('{}/histogram_group3_prior&posterior.png'.format(path))

def histogramspriorposterior_logscale_top_kl_div(posteriorsamples, prior_std_dev, path, color):
    # histogram with 95th percentile of prior
    prior_estimatedvariables = np.load('posteriorhistograms/prior_estimatedvariables.npy')
    variable_dictionary = np.load('posteriorhistograms/dictionary_estimatedvariables.npy')
    std_prior = prior_std_dev
    k24_index = np.where(variable_dictionary == 'k24')[0][0]

    plt.figure(figsize=(5, 9))
    for i in np.arange(0, 33):
        # plt.subplots_adjust(bottom=0.05, top=0.99, left=0.05, right=0.99)
        parameter = np.log(posteriorsamples[:, i])
        mu = np.log(prior_estimatedvariables[i])
        prior = np.log(np.random.lognormal(mean=mu, sigma=std_prior, size=len(parameter)))
        if variable_dictionary[i] == 'deg_ratio' or variable_dictionary[i] == 'k27' or variable_dictionary[i] == 'k24':
            if variable_dictionary[i] == 'deg_ratio':
                j = 2
            if variable_dictionary[i] == 'k27':
                j = 1
            if variable_dictionary[i] == 'k24':
                j = 3
            plt.subplot(3, 1, j)
            plt.subplots_adjust(bottom=0.1, top=0.99, left=0.05, right=0.95)
            sns.histplot(data=parameter, color='{}'.format(color))
            sns.histplot(data=prior, color='{}'.format('slategrey'), kde=True)
            if variable_dictionary[i] == 'deg_ratio':
                plt.legend(['prior','posterior'], fontsize=14)
            plt.xlabel('ln({})'.format(variable_dictionary[i]), fontsize=14)
            plt.yticks([0, 25000], '')
            plt.ylabel('')
            plt.ylim([0, 25000])
            #plt.tight_layout(pad=0.2, h_pad=0.02, w_pad=0.02)
            plt.savefig('{}/histogram_top_kl_div_prior&posterior.png'.format(path), dpi=300)
