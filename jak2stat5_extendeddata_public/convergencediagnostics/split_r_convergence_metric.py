import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import statsmodels.api as sm


def qqplot(posteriorsamples):
    variable_dictionary = np.load('convergencediagnostics/dictionary_estimatedvariables.npy')

    fig = plt.figure(figsize=(13, 10))
    for i in np.arange(0, 11):
        parameter = np.subtract(np.log(posteriorsamples[:, i]), np.mean(np.log(posteriorsamples[:, i])))
        ax = fig.add_subplot(3, 4, i + 1)
        plt.subplots_adjust(bottom=0.1, top=0.99, left=0.05, right=0.95)
        sm.qqplot(parameter, fit=True, line='45', ax=ax, markerfacecolor = 'slategrey', c = 'slategrey')
        plt.legend(['ln({})'.format(variable_dictionary[i])])
        plt.xlabel('')
        plt.ylabel('')
        if i == 0:
            plt.ylabel('Sample Quantiles')
            plt.xlabel('Theoretical Quantiles')
        plt.xlim([-3.1, 3.1])
        plt.xticks(np.arange(-3,4))
        plt.ylim([-3.1, 3.1])
        plt.yticks(np.arange(-3, 4))
    plt.savefig('convergencediagnostics/qq_plots_ln(posterior)/Q_Q_Plot_group1.png', dpi=150)
    plt.close(fig)

    fig = plt.figure(figsize=(13, 10))
    for i in np.arange(11, 22):
        parameter = np.subtract(np.log(posteriorsamples[:, i]), np.mean(np.log(posteriorsamples[:, i])))
        ax = fig.add_subplot(3, 4, i - 11 + 1)
        plt.subplots_adjust(bottom=0.1, top=0.99, left=0.05, right=0.95)
        sm.qqplot(parameter, fit=True, line='45', ax=ax, markerfacecolor = 'slategrey', c = 'slategrey')
        plt.legend(['ln({})'.format(variable_dictionary[i])])
        plt.xlabel('Sample Quantiles')
        plt.ylabel('Theoretical Quantiles')
        plt.xlabel('')
        plt.ylabel('')
        if i == 0+11:
            plt.ylabel('Sample Quantiles')
            plt.xlabel('Theoretical Quantiles')
        plt.xlim([-3.1, 3.1])
        plt.xticks(np.arange(-3,4))
        plt.ylim([-3.1, 3.1])
        plt.yticks(np.arange(-3, 4))
    plt.savefig('convergencediagnostics/qq_plots_ln(posterior)/Q_Q_Plot_group2.png', dpi=150)
    plt.close(fig)

    fig = plt.figure(figsize=(13, 10))
    for i in np.arange(22, 33):
        parameter = np.subtract(np.log(posteriorsamples[:, i]), np.mean(np.log(posteriorsamples[:, i])))
        ax = fig.add_subplot(3, 4, i - 22 + 1)
        plt.subplots_adjust(bottom=0.1, top=0.99, left=0.05, right=0.95)
        sm.qqplot(parameter, fit=True, line='45', ax=ax, markerfacecolor = 'slategrey', c = 'slategrey')
        plt.legend(['ln({})'.format(variable_dictionary[i])])
        plt.xlabel('Sample Quantiles')
        plt.ylabel('Theoretical Quantiles')
        plt.xlabel('')
        plt.ylabel('')
        if i == 0+22:
            plt.ylabel('Sample Quantiles')
            plt.xlabel('Theoretical Quantiles')
        plt.xlim([-3.1, 3.1])
        plt.xticks(np.arange(-3,4))
        plt.ylim([-3.1, 3.1])
        plt.yticks(np.arange(-3, 4))
    plt.savefig('convergencediagnostics/qq_plots_ln(posterior)/Q_Q_Plot_group3.png', dpi=150)
    plt.close(fig)

def split_r_convergence_metric(estimated_parameters1, estimated_parameters2, estimated_parameters3):
    gelman_rubin = []
    variable_dictionary = np.load('convergencediagnostics/dictionary_estimatedvariables.npy')
    for i in np.arange(0, 33):
        N = int(len(estimated_parameters1[:, i])/2)
        M = 6
        'Gelman Rubin'
        chain1 = np.log(estimated_parameters1[0:int(len(estimated_parameters1[:, i])/2), i])
        chain2 = np.log(estimated_parameters1[int(len(estimated_parameters1[:, i])/2):-1, i])
        chain3 = np.log(estimated_parameters2[0:int(len(estimated_parameters2[:, i])/2), i])
        chain4 = np.log(estimated_parameters2[int(len(estimated_parameters2[:, i])/2):-1, i])
        chain5 = np.log(estimated_parameters3[0:int(len(estimated_parameters3[:, i])/2), i])
        chain6 = np.log(estimated_parameters3[int(len(estimated_parameters3[:, i])/2):-1, i])
        mean1 = np.mean(chain1)
        mean2 = np.mean(chain2)
        mean3 = np.mean(chain3)
        mean4 = np.mean(chain4)
        mean5 = np.mean(chain5)
        mean6 = np.mean(chain6)
        averagemean = np.mean(np.array([mean1, mean2, mean3, mean4, mean5, mean6]))
        var1 = np.var(chain1)
        var2 = np.var(chain2)
        var3 = np.var(chain3)
        var4 = np.var(chain4)
        var5 = np.var(chain5)
        var6 = np.var(chain6)
        B_var = []
        for k in [mean1, mean2, mean3, mean4, mean5, mean6]:
            B_vark = (k - averagemean) ** 2
            B_var.append(B_vark)
        B_var = np.array(B_var)
        B = (N / (M - 1)) * np.sum(B_var)
        W = (1 / M) * np.sum(np.array([var1, var2, var3, var4, var5, var6]))
        posterior_variance = ((N - 1) / N) * W + (1 / N) * B
        R = np.around(np.sqrt(posterior_variance / W), 2)
        gelman_rubin.append(R)
        print(variable_dictionary[i], R)

    # Gelman Rubin Histogram
    fig = plt.figure()
    sns.histplot(data=gelman_rubin, color='dimgrey', alpha=0.75)
    plt.axvline(x=1.1, ymin=0, ymax=1, label='convergence cutoff', color='firebrick', linewidth=2)
    plt.xlabel('Split-Potential Scale Reduction Factor')
    plt.xticks([0.975, 1, 1.025, 1.050, 1.075, 1.1])
    plt.savefig('convergencediagnostics/split_gelman_rubin_histogram.png', dpi=300)