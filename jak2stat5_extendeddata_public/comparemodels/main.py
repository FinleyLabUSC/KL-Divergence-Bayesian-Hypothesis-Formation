import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

#run script to generate figure 6 and supplementary figure S10

fontsize = 14
# stacked bar plot, 50% cutoff
baseline_sustained = 13+71
baseline_transient = 11+5
dataint_sustained = 6+4
dataint_transient = 51+38

color = ['silver', 'slategrey']
labels = ['Original \n data', 'Extended \n data']

p1 = plt.bar(labels, np.array([baseline_transient, dataint_transient]), color=['silver','lightcoral'])
p2 = plt.bar(labels, np.array([baseline_sustained, dataint_sustained]), color=['dimgrey','firebrick'], bottom=np.array([baseline_transient, dataint_transient]))
plt.ylabel('Long-Term SOCS Behavior Frequency', fontsize=fontsize)
plt.xticks([0, 1], labels, fontsize=fontsize)
plt.show()
plt.savefig('comparemodels/stacked_bar_long_term_frequency_50%.png', dpi=300)

# stacked bar plot, 30% cutoff
baseline_sustained = 16+78
baseline_transient = 4+1
dataint_sustained = 17+20
dataint_transient = 38+25

color = ['silver', 'slategrey']
labels = ['Original \n data', 'Extended \n data']

p1 = plt.bar(labels, np.array([baseline_transient, dataint_transient]), color=['silver','lightcoral'])
p2 = plt.bar(labels, np.array([baseline_sustained, dataint_sustained]), color=['dimgrey','firebrick'], bottom=np.array([baseline_transient, dataint_transient]))
plt.ylabel('Long-Term SOCS Behavior Frequency', fontsize=fontsize)
plt.xticks([0, 1], labels, fontsize=fontsize)
plt.show()
plt.savefig('comparemodels/stacked_bar_long_term_frequency_30%.png', dpi=300)

# stacked bar plot, 70% cutoff
baseline_sustained = 8+58
baseline_transient = 24+9
dataint_sustained = 1+1
dataint_transient = 56+41

color = ['silver', 'slategrey']
labels = ['Original \n data', 'Extended \n data']

p1 = plt.bar(labels, np.array([baseline_transient, dataint_transient]), color=['silver','lightcoral'])
p2 = plt.bar(labels, np.array([baseline_sustained, dataint_sustained]), color=['dimgrey','firebrick'], bottom=np.array([baseline_transient, dataint_transient]))
plt.ylabel('Long-Term SOCS Behavior Frequency', fontsize=fontsize)
plt.xticks([0, 1], labels, fontsize=fontsize)
plt.show()
plt.savefig('comparemodels/stacked_bar_long_term_frequency_70%.png', dpi=300)

#k24 posterior comparison
posterior_samples_dataint = np.load('generateinputs/all_chains_baseline and timepoint_proposal 0.1_thinnedestimatedvariables.npy')
posterior_samples_baseline = np.load('comparemodels/all_chains_baseline_proposal 0.1_thinnedestimatedvariables.npy')
variable_dictionary = np.load('comparemodels/dictionary_estimatedvariables.npy')
k24_index = np.where(variable_dictionary == 'k24')[0][0]
dataint_k24_samples = np.log(posterior_samples_dataint[:, k24_index])
baseline_k24_samples = np.log(posterior_samples_baseline[:, k24_index])

sns.histplot(baseline_k24_samples, color='black')
sns.histplot(dataint_k24_samples, color='firebrick')
plt.xlabel('ln(SOCS:RECEPTOR Degradation Rate)', fontsize=fontsize)
plt.ylabel('Count', fontsize=fontsize)
#plt.legend(['baseline data', 'extended data'])
plt.savefig('comparemodels/histogram_posterior_comparison.png', dpi=300)

# comparison of all refined posteriors
# histogram with 95th percentile of prior
prior_estimatedvariables = np.load('posteriorhistograms/prior_estimatedvariables.npy')
variable_dictionary = np.load('posteriorhistograms/dictionary_estimatedvariables.npy')
k24_index = np.where(variable_dictionary == 'k24')[0][0]
std_prior = np.sqrt(2)
path = 'comparemodels'

plt.figure(figsize=(13, 10))
for i in np.arange(0, 11):
    parameter = np.log(posterior_samples_dataint[:, i])
    parameter_original = np.log(posterior_samples_baseline[:, i])
    mu = np.log(prior_estimatedvariables[i])
    #prior = np.log(np.random.lognormal(mean=mu, sigma=std_prior, size=len(parameter)))
    upper_percentile = mu + 2 * std_prior
    lower_percentile = mu - 2 * std_prior
    plt.subplot(3, 4, i + 1)
    plt.subplots_adjust(bottom=0.05, top=0.95, left=0.05, right=0.99)
    sns.histplot(data=parameter_original, color='black')
    sns.histplot(data=parameter, color='firebrick')
    #sns.histplot(data=prior, color='slategrey', kde=True)
    plt.vlines(x=upper_percentile, ymin=0, ymax=1000, color='slategrey', linestyles='dashed')
    plt.vlines(x=lower_percentile, ymin=0, ymax=1000, color='slategrey', linestyles='dashed')
    plt.xlabel('ln({})'.format(variable_dictionary[i]), fontsize=14)
    if i != 0 and i != 4 and i != 8:
        plt.yticks([0, 1000], '')
    plt.ylabel('', fontsize=14)
    plt.ylim([0, 1000])
plt.savefig('{}/histogram_group1_posterior&posterior.png'.format(path))

plt.figure(figsize=(13, 10))
for i in np.arange(11, 22):
    parameter = np.log(posterior_samples_dataint[:, i])
    parameter_original = np.log(posterior_samples_baseline[:, i])
    mu = np.log(prior_estimatedvariables[i])
    #prior = np.log(np.random.lognormal(mean=mu, sigma=std_prior, size=len(parameter)))
    upper_percentile = mu + 2 * std_prior
    lower_percentile = mu - 2 * std_prior
    plt.subplot(3, 4, i - 11 + 1)
    plt.subplots_adjust(bottom=0.05, top=0.95, left=0.05, right=0.99)
    if i == k24_index:
        sns.histplot(data=parameter_original, color='black')
        sns.histplot(data=parameter, color='firebrick')
        #sns.histplot(data=prior, color='black', kde=True)
    else:
        sns.histplot(data=parameter_original, color='black')
        sns.histplot(data=parameter, color='firebrick')
        #sns.histplot(data=prior, color='slategrey', kde=True)
    plt.vlines(x=upper_percentile, ymin=0, ymax=1000, color='slategrey', linestyles='dashed')
    plt.vlines(x=lower_percentile, ymin=0, ymax=1000, color='slategrey', linestyles='dashed')
    plt.xlabel('ln({})'.format(variable_dictionary[i]), fontsize=14)
    if i != 0+11 and i != 4+11 and i != 8+11:
        plt.yticks([0, 1000], '')
    plt.ylabel('', fontsize=14)
    plt.ylim([0, 1000])
plt.savefig('{}/histogram_group2_prior&posterior.png'.format(path))

plt.figure(figsize=(13, 10))
for i in np.arange(22, 33):
    parameter = np.log(posterior_samples_dataint[:, i])
    parameter_original = np.log(posterior_samples_baseline[:, i])
    mu = np.log(prior_estimatedvariables[i])
    #prior = np.log(np.random.lognormal(mean=mu, sigma=std_prior, size=len(parameter)))
    upper_percentile = mu + 2 * std_prior
    lower_percentile = mu - 2 * std_prior
    plt.subplot(3, 4, i - 22 + 1)
    plt.subplots_adjust(bottom=0.05, top=0.95, left=0.05, right=0.99)
    sns.histplot(data=parameter_original, color='black')
    sns.histplot(data=parameter, color='firebrick')
    #sns.histplot(data=prior, color='slategrey', kde=True)
    plt.vlines(x=upper_percentile, ymin=0, ymax=1000, color='slategrey', linestyles='dashed')
    plt.vlines(x=lower_percentile, ymin=0, ymax=1000, color='slategrey', linestyles='dashed')
    plt.xlabel('ln({})'.format(variable_dictionary[i]), fontsize=14)
    if i != 0 and i != 4 and i != 8:
        plt.yticks([0, 1000], '')
    plt.ylabel('', fontsize=14)
    plt.ylim([0, 1000])
plt.savefig('{}/histogram_group3_prior&posterior.png'.format(path))