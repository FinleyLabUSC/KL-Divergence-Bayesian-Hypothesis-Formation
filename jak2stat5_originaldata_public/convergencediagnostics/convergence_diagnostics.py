import numpy as np
from pandas import Series
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

def convergencediagonistics(estimated_parameters1, estimated_parameters2, estimated_parameters3, estimated_parameters4, estimated_parameters5 ):
    chaincount = 5
    gelman_rubin = []
    variable_dictionary = np.load('convergencediagnostics/dictionary_estimatedvariables.npy')
    estimated_parametersinit = np.load('convergencediagnostics/prior_estimatedvariables.npy')
    for i in np.arange(0, 2):
        N = int(len(estimated_parameters1[:, i])/2)
        M = 10
        'Gelman Rubin'
        chain1 = np.log(estimated_parameters1[0:int(len(estimated_parameters1[:, i])/2), i])
        chain2 = np.log(estimated_parameters1[int(len(estimated_parameters1[:, i])/2):-1, i])
        chain3 = np.log(estimated_parameters2[0:int(len(estimated_parameters2[:, i])/2), i])
        chain4 = np.log(estimated_parameters2[int(len(estimated_parameters2[:, i])/2):-1, i])
        chain5 = np.log(estimated_parameters3[0:int(len(estimated_parameters3[:, i])/2), i])
        chain6 = np.log(estimated_parameters3[int(len(estimated_parameters3[:, i])/2):-1, i])
        chain7 = np.log(estimated_parameters4[0:int(len(estimated_parameters4[:, i])/2), i])
        chain8 = np.log(estimated_parameters4[int(len(estimated_parameters4[:, i])/2):-1, i])
        chain9 = np.log(estimated_parameters5[0:int(len(estimated_parameters5[:, i])/2), i])
        chain10 = np.log(estimated_parameters5[int(len(estimated_parameters5[:, i])/2):-1, i])
        mean1 = np.mean(chain1)
        mean2 = np.mean(chain2)
        mean3 = np.mean(chain3)
        mean4 = np.mean(chain4)
        mean5 = np.mean(chain5)
        mean6 = np.mean(chain6)
        mean7 = np.mean(chain7)
        mean8 = np.mean(chain8)
        mean9 = np.mean(chain9)
        mean10 = np.mean(chain10)
        averagemean = np.mean(np.array([mean1, mean2, mean3, mean4, mean5, mean6, mean7, mean8, mean9, mean10]))
        var1 = np.var(chain1)
        var2 = np.var(chain2)
        var3 = np.var(chain3)
        var4 = np.var(chain4)
        var5 = np.var(chain5)
        var6 = np.var(chain6)
        var7 = np.var(chain3)
        var8 = np.var(chain4)
        var9 = np.var(chain5)
        var10 = np.var(chain6)
        B_var = []
        for k in [mean1, mean2, mean3, mean4, mean5, mean6, mean7, mean8, mean9, mean10]:
            B_vark = (k - averagemean) ** 2
            B_var.append(B_vark)
        B_var = np.array(B_var)
        B = (N / (M - 1)) * np.sum(B_var)
        W = (1 / M) * np.sum(np.array([var1, var2, var3, var4, var5, var6, var7, var8, var9, var10]))
        posterior_variance = ((N - 1) / N) * W + (1 / N) * B
        R = np.around(np.sqrt(posterior_variance / W), 2)
        gelman_rubin.append(R)
        print(variable_dictionary[i], R)

        'Trace Plot'
        fig = plt.figure(figsize=(10, 7))
        plt.subplots_adjust(bottom=0.01, top=0.95)
        plt.subplot(4, 1, 1)
        plt.plot(np.arange(0, len(estimated_parameters1)), estimated_parameters1[:, i], color='blue', linewidth=0.5)
        plt.plot(np.arange(0, len(estimated_parameters2)), estimated_parameters2[:, i], color='red', linewidth=0.5)
        plt.plot(np.arange(0, len(estimated_parameters3)), estimated_parameters3[:, i], color='orange', linewidth=0.5)
        plt.plot(np.arange(0, len(estimated_parameters4)), estimated_parameters4[:, i], color='purple', linewidth=0.5)
        plt.plot(np.arange(0, len(estimated_parameters5)), estimated_parameters5[:, i], color='green', linewidth=0.5)
        plt.title('Convergence Diagnostics for {}: Trace, Histogram, Autocorrelation, Split-R: {}'.format(variable_dictionary[i], R))
        plt.ylabel('Value', fontsize=8)
        plt.ylim([-estimated_parametersinit[i], estimated_parametersinit[i] * 1e2])
        plt.yticks(fontsize=8)
        plt.xticks(fontsize=8)
        plt.legend(list(np.arange(0, chaincount)+1), loc='upper right')

        'Histogram'
        plt.subplot(4, 1, 2)
        sns.histplot(data=estimated_parameters1[:, i], color='blue')
        sns.histplot(data=estimated_parameters2[:, i], color='red')
        sns.histplot(data=estimated_parameters3[:, i], color='orange')
        sns.histplot(data=estimated_parameters4[:, i], color='purple')
        sns.histplot(data=estimated_parameters5[:, i], color='green')
        plt.legend(list(np.arange(0, chaincount)+1), loc='upper right')
        plt.xlim([-estimated_parametersinit[i], estimated_parametersinit[i] * 1e2])
        plt.ylim([0, 20000])
        plt.yticks(fontsize=8)
        plt.xticks(fontsize=8)
        plt.ylabel('Count', fontsize=8)

        'Autocorrelation'
        plt.subplot(4, 1, 3)
        pearson1 = []
        maxlag = 20000
        interval = int(maxlag / 1000)
        for j in np.arange(0, maxlag, interval):
            autocorrelation1 = Series.autocorr(pd.Series(estimated_parameters1[:, i]), lag=j)
            pearson1.append(autocorrelation1)
        pearson2 = []
        for j in np.arange(0, maxlag, interval):
            autocorrelation2 = Series.autocorr(pd.Series(estimated_parameters2[:, i]), lag=j)
            pearson2.append(autocorrelation2)
        pearson3 = []
        for j in np.arange(0, maxlag, interval):
            autocorrelation3 = Series.autocorr(pd.Series(estimated_parameters3[:, i]), lag=j)
            pearson3.append(autocorrelation3)
        pearson4 = []
        for j in np.arange(0, maxlag, interval):
            autocorrelation4 = Series.autocorr(pd.Series(estimated_parameters4[:, i]), lag=j)
            pearson4.append(autocorrelation4)
        pearson5 = []
        for j in np.arange(0, maxlag, interval):
            autocorrelation5 = Series.autocorr(pd.Series(estimated_parameters5[:, i]), lag=j)
            pearson5.append(autocorrelation5)

        plt.plot(np.arange(0, maxlag, interval), pearson1, color='blue')
        plt.fill_between(np.arange(0, maxlag, interval), pearson1, 0, color='blue')
        plt.plot(np.arange(0, maxlag, interval), pearson2, color='red', alpha=0.5)
        plt.fill_between(np.arange(0, maxlag, interval), pearson2, 0, color='red', alpha=0.5)
        plt.plot(np.arange(0, maxlag, interval), pearson3, color='orange', alpha=0.5)
        plt.fill_between(np.arange(0, maxlag, interval), pearson3, 0, color='orange', alpha=0.5)
        plt.plot(np.arange(0, maxlag, interval), pearson4, color='purple')
        plt.fill_between(np.arange(0, maxlag, interval), pearson4, 0, color='purple')
        plt.plot(np.arange(0, maxlag, interval), pearson5, color='green')
        plt.fill_between(np.arange(0, maxlag, interval), pearson5, 0, color='green')
        plt.legend(list(np.arange(0, chaincount)+1), loc='upper right')
        plt.ylabel('Autocorrelation', fontsize=8)
        plt.xlabel('Lag', fontsize=8)
        plt.ylim([-1.25, 1.25])
        plt.yticks(fontsize=8)
        plt.xticks(fontsize=8)

        'Chain Stats'
        ax = plt.subplot(4, 1, 4)
        chain1 = estimated_parameters1[:, i]
        chain2 = estimated_parameters2[:, i]
        chain3 = estimated_parameters3[:, i]
        chain4 = estimated_parameters4[:, i]
        chain5 = estimated_parameters5[:, i]
        chains = np.transpose(np.array([chain1, chain2[0:len(chain1)], chain3[0:len(chain1)], chain4, chain5]))
        df = pd.DataFrame(chains, columns=list(np.arange(0, chaincount)+1))
        describe = df.describe()
        table = ax.table(cellText=describe.values, colLabels=describe.columns, rowLabels=describe.index, loc='center',
                         fontsize='x-small')
        table.set_fontsize(8)
        ax.axis('off')
        ax.axis('tight')

        'Save'
        plt.savefig('convergencediagnostics/mcmc_qualitative_diagnostics/diagnostics_for_{}.png'.format(variable_dictionary[i]))
        plt.close(fig)

def convergencediagonistics_logscale(estimated_parameters1, estimated_parameters2, estimated_parameters3, estimated_parameters4, estimated_parameters5 ):
    fontsize=12
    chaincount = 5
    gelman_rubin = []
    variable_dictionary = np.load('convergencediagnostics/dictionary_estimatedvariables.npy')
    estimated_parametersinit = np.load('convergencediagnostics/prior_estimatedvariables.npy')
    for i in np.arange(0, 33):
        N = int(len(estimated_parameters1[:, i])/2)
        M = 10
        'Gelman Rubin'
        chain1 = np.log(estimated_parameters1[0:int(len(estimated_parameters1[:, i])/2), i])
        chain2 = np.log(estimated_parameters1[int(len(estimated_parameters1[:, i])/2):-1, i])
        chain3 = np.log(estimated_parameters2[0:int(len(estimated_parameters2[:, i])/2), i])
        chain4 = np.log(estimated_parameters2[int(len(estimated_parameters2[:, i])/2):-1, i])
        chain5 = np.log(estimated_parameters3[0:int(len(estimated_parameters3[:, i])/2), i])
        chain6 = np.log(estimated_parameters3[int(len(estimated_parameters3[:, i])/2):-1, i])
        chain7 = np.log(estimated_parameters4[0:int(len(estimated_parameters4[:, i])/2), i])
        chain8 = np.log(estimated_parameters4[int(len(estimated_parameters4[:, i])/2):-1, i])
        chain9 = np.log(estimated_parameters5[0:int(len(estimated_parameters5[:, i])/2), i])
        chain10 = np.log(estimated_parameters5[int(len(estimated_parameters5[:, i])/2):-1, i])
        mean1 = np.mean(chain1)
        mean2 = np.mean(chain2)
        mean3 = np.mean(chain3)
        mean4 = np.mean(chain4)
        mean5 = np.mean(chain5)
        mean6 = np.mean(chain6)
        mean7 = np.mean(chain7)
        mean8 = np.mean(chain8)
        mean9 = np.mean(chain9)
        mean10 = np.mean(chain10)
        averagemean = np.mean(np.array([mean1, mean2, mean3, mean4, mean5, mean6, mean7, mean8, mean9, mean10]))
        var1 = np.var(chain1)
        var2 = np.var(chain2)
        var3 = np.var(chain3)
        var4 = np.var(chain4)
        var5 = np.var(chain5)
        var6 = np.var(chain6)
        var7 = np.var(chain3)
        var8 = np.var(chain4)
        var9 = np.var(chain5)
        var10 = np.var(chain6)
        B_var = []
        for k in [mean1, mean2, mean3, mean4, mean5, mean6, mean7, mean8, mean9, mean10]:
            B_vark = (k - averagemean) ** 2
            B_var.append(B_vark)
        B_var = np.array(B_var)
        B = (N / (M - 1)) * np.sum(B_var)
        W = (1 / M) * np.sum(np.array([var1, var2, var3, var4, var5, var6, var7, var8, var9, var10]))
        posterior_variance = ((N - 1) / N) * W + (1 / N) * B
        R = np.around(np.sqrt(posterior_variance / W), 2)
        gelman_rubin.append(R)
        print(variable_dictionary[i], R)

        'Trace Plot'
        fig = plt.figure(figsize=(10, 7))
        plt.subplots_adjust(bottom=0.01, top=0.95)
        plt.subplot(3, 1, 1)
        plt.plot(np.arange(0, len(estimated_parameters1)), np.log(estimated_parameters1[:, i]), color='blue', linewidth=0.5)
        plt.plot(np.arange(0, len(estimated_parameters2)), np.log(estimated_parameters2[:, i]), color='red', linewidth=0.5)
        plt.plot(np.arange(0, len(estimated_parameters3)), np.log(estimated_parameters3[:, i]), color='orange', linewidth=0.5)
        plt.plot(np.arange(0, len(estimated_parameters4)), np.log(estimated_parameters4[:, i]), color='purple', linewidth=0.5)
        plt.plot(np.arange(0, len(estimated_parameters5)), np.log(estimated_parameters5[:, i]), color='green', linewidth=0.5)
        plt.title('Convergence Diagnostics for {}: Trace, Histogram, Autocorrelation, Split-R: {}'.format(variable_dictionary[i], R))
        plt.ylabel('ln({})'.format(variable_dictionary[i]), fontsize=8)
        #plt.ylim([np.log(estimated_parametersinit[i]), np.log(estimated_parametersinit[i] * 1e2)])
        plt.yticks(fontsize=fontsize)
        plt.xticks(fontsize=fontsize)
        plt.legend(list(np.arange(0, chaincount)+1), loc='upper right')
        plt.xlabel('chain sample', fontsize=fontsize)

        'Histogram'
        plt.subplot(3, 1, 2)
        sns.histplot(data=np.log(estimated_parameters1[:, i]), color='blue', element="step", fill=False)
        sns.histplot(data=np.log(estimated_parameters2[:, i]), color='red', element="step", fill=False)
        sns.histplot(data=np.log(estimated_parameters3[:, i]), color='orange', element="step", fill=False)
        sns.histplot(data=np.log(estimated_parameters4[:, i]), color='purple', element="step", fill=False)
        sns.histplot(data=np.log(estimated_parameters5[:, i]), color='green', element="step", fill=False)
        plt.legend(list(np.arange(0, chaincount)+1), loc='upper right')
        #plt.xlim([-estimated_parametersinit[i], estimated_parametersinit[i] * 1e2])
        plt.ylim([0, 10000])
        plt.yticks(fontsize=fontsize)
        plt.xticks(fontsize=fontsize)
        plt.ylabel('Count', fontsize=fontsize)
        plt.xlabel('ln({})'.format(variable_dictionary[i]), fontsize=fontsize)

        'Autocorrelation'
        plt.subplot(3, 1, 3)
        pearson1 = []
        maxlag = 20000
        interval = int(maxlag / 1000)
        for j in np.arange(0, maxlag, interval):
            autocorrelation1 = Series.autocorr(pd.Series(estimated_parameters1[:, i]), lag=j)
            pearson1.append(autocorrelation1)
        pearson2 = []
        for j in np.arange(0, maxlag, interval):
            autocorrelation2 = Series.autocorr(pd.Series(estimated_parameters2[:, i]), lag=j)
            pearson2.append(autocorrelation2)
        pearson3 = []
        for j in np.arange(0, maxlag, interval):
            autocorrelation3 = Series.autocorr(pd.Series(estimated_parameters3[:, i]), lag=j)
            pearson3.append(autocorrelation3)
        pearson4 = []
        for j in np.arange(0, maxlag, interval):
            autocorrelation4 = Series.autocorr(pd.Series(estimated_parameters4[:, i]), lag=j)
            pearson4.append(autocorrelation4)
        pearson5 = []
        for j in np.arange(0, maxlag, interval):
            autocorrelation5 = Series.autocorr(pd.Series(estimated_parameters5[:, i]), lag=j)
            pearson5.append(autocorrelation5)

        plt.plot(np.arange(0, maxlag, interval), pearson1, color='blue')
        plt.fill_between(np.arange(0, maxlag, interval), pearson1, 0, color='blue', alpha=0.3)
        plt.plot(np.arange(0, maxlag, interval), pearson2, color='red', alpha=0.5)
        plt.fill_between(np.arange(0, maxlag, interval), pearson2, 0, color='red', alpha=0.3)
        plt.plot(np.arange(0, maxlag, interval), pearson3, color='orange', alpha=0.5)
        plt.fill_between(np.arange(0, maxlag, interval), pearson3, 0, color='orange', alpha=0.3)
        plt.plot(np.arange(0, maxlag, interval), pearson4, color='purple')
        plt.fill_between(np.arange(0, maxlag, interval), pearson4, 0, color='purple', alpha=0.3)
        plt.plot(np.arange(0, maxlag, interval), pearson5, color='green')
        plt.fill_between(np.arange(0, maxlag, interval), pearson5, 0, color='green', alpha=0.3)
        plt.legend(list(np.arange(0, chaincount)+1), loc='upper right')
        plt.ylabel('Autocorrelation', fontsize=fontsize)
        plt.xlabel('Lag', fontsize=fontsize)
        plt.ylim([-1.25, 1.25])
        plt.yticks(fontsize=fontsize)
        plt.xticks(fontsize=fontsize)

        plt.tight_layout()

        'Save'
        plt.savefig('convergencediagnostics/mcmc_qualitative_diagnostics/diagnostics_for_{}.png'.format(variable_dictionary[i]), dpi=300)
        plt.close(fig)
