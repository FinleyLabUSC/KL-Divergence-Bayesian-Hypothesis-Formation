import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import uniform
from kde import estimate_univariate_kde
from kl_divergence import estimate_kl_divergence_posterior_loguninform_prior_totsch

#load data
posteriors_DSAS = pd.read_hdf('ER_noAD_iid.h5')
posteriors_DREG = pd.read_hdf('reg_noAD_iid.h5')

DSAS_parameters = posteriors_DSAS.columns
DREG_parameters = posteriors_DREG.columns

# kdes for DSAS
estimate_univariate_kde(posteriorsamples=np.array(posteriors_DSAS), name='DSAS')

# kdes for DREG
estimate_univariate_kde(posteriorsamples=np.array(posteriors_DREG), name='DREG')

#kld for DSAS kd
KDE_DSAS = np.load('kdes_DSAS.npy', allow_pickle=True)
KDE_kd = KDE_DSAS[np.argwhere(np.array(DSAS_parameters)=='kd')][0][0]
estimate_kl_divergence_posterior_loguninform_prior_totsch(posteriorsamples_all=np.array(posteriors_DSAS["kd"]), prior_lower=-8, prior_upper=1, KDE=KDE_kd, name='kd_DSAS')
kd_KLD = np.load('kullback_leibler_divergence_kd_DSAS.npy')

#kld for DSAS k_reg/k_BR
parameter = 'k_reg'
KDE_DSAS = np.load('kdes_DSAS.npy', allow_pickle=True)
KDE_kd = KDE_DSAS[np.argwhere(np.array(DSAS_parameters)=="{}".format(parameter))][0][0]
estimate_kl_divergence_posterior_loguninform_prior_totsch(posteriorsamples_all=np.array(posteriors_DSAS["{}".format(parameter)]), prior_lower=-0.5, prior_upper=9, KDE=KDE_kd, name="{}_DSAS".format(parameter))
k_BR_KLD_DSAS = np.load('kullback_leibler_divergence_{}_DSAS.npy'.format(parameter))

#kld for DSAS S0
parameter = 'S0'
KDE_DSAS = np.load('kdes_DSAS.npy', allow_pickle=True)
KDE_kd = KDE_DSAS[np.argwhere(np.array(DSAS_parameters)=="{}".format(parameter))][0][0]
estimate_kl_divergence_posterior_loguninform_prior_totsch(posteriorsamples_all=np.array(posteriors_DSAS["{}".format(parameter)]), prior_lower=-12, prior_upper=-3, KDE=KDE_kd, name="{}_DSAS".format(parameter))
S0_KLD_DSAS = np.load('kullback_leibler_divergence_{}_DSAS.npy'.format(parameter))

#kld for DSAS KA
parameter = 'KA'
KDE_DSAS = np.load('kdes_DSAS.npy', allow_pickle=True)
KDE_kd = KDE_DSAS[np.argwhere(np.array(DSAS_parameters)=="{}".format(parameter))][0][0]
estimate_kl_divergence_posterior_loguninform_prior_totsch(posteriorsamples_all=np.array(posteriors_DSAS["{}".format(parameter)]), prior_lower=3, prior_upper=12, KDE=KDE_kd, name="{}_DSAS".format(parameter))
KA_KLD_DSAS = np.load('kullback_leibler_divergence_{}_DSAS.npy'.format(parameter))

#kld for DSAS_k'ER/k
parameter = 'k'
KDE_DSAS = np.load('kdes_DSAS.npy', allow_pickle=True)
KDE_kd = KDE_DSAS[np.argwhere(np.array(DSAS_parameters)=="{}".format(parameter))][0][0]
estimate_kl_divergence_posterior_loguninform_prior_totsch(posteriorsamples_all=np.array(posteriors_DSAS["{}".format(parameter)]), prior_lower=-5.15, prior_upper=-0.15, KDE=KDE_kd, name="{}_DSAS".format(parameter))
k_KLD_DSAS = np.load('kullback_leibler_divergence_{}_DSAS.npy'.format(parameter))

#KLD plot
kd_KLD = np.load('kullback_leibler_divergence_kd.npy')
k_BR_KLD_DSAS = np.load('kullback_leibler_divergence_k_reg_DSAS.npy')
S0_KLD_DSAS = np.load('kullback_leibler_divergence_S0_DSAS.npy')
KA_KLD_DSAS = np.load('kullback_leibler_divergence_KA_DSAS.npy')
k_KLD_DSAS = np.load('kullback_leibler_divergence_k_DSAS.npy')
names = ["kd", "KA", "S0", "k'_ER", "k_BR"]
kld = [kd_KLD, KA_KLD_DSAS, S0_KLD_DSAS, k_KLD_DSAS, k_BR_KLD_DSAS]
plt.scatter(np.arange(1, 6), kld, color='black')
plt.xlabel("Parameters", fontsize=14)
plt.xticks(np.arange(1, 6), names)
plt.ylabel("KL Divergence", fontsize=14)
plt.savefig('KLD_Totsch.png', dpi=300)
plt.close()

#histogram plot
plt.figure(figsize=(5, 9))
plt.subplot(3, 1, 1)
sns.kdeplot(np.array(posteriors_DSAS["k_reg"]), color='black', fill=True, alpha=1)
x = np.linspace(-0.51, 9.01, 1000)
plt.plot(x, uniform.pdf(x, loc=-0.5, scale=(9--0.5)), color="slategrey")
plt.fill_between(x=[-0.5, 9], y1=uniform.pdf(x[3], loc=-0.5, scale=(9--0.5)), y2=0, color="slategrey", alpha=0.6)
plt.xlabel('log(k_BR)', fontsize=14)
plt.yticks([0, 1], labels=['',''])
plt.ylabel('')
plt.subplot(3, 1, 2)
sns.kdeplot(np.array(posteriors_DSAS["k"]), color='black', fill=True, alpha=1, label='posterior')
x = np.linspace(-5.16, -0.14, 1000)
plt.plot(x, uniform.pdf(x, loc=-5.15, scale=(-0.15--5.15)), color="slategrey", label='prior')
plt.fill_between(x=[-5.15, -0.15], y1=uniform.pdf(x[3], loc=-5.15, scale=(-0.15--5.15)), y2=0, color="slategrey", alpha=0.6)
plt.xlabel('log(k_ER)', fontsize=14)
plt.yticks([0, 1], labels=['',''])
plt.ylabel('')
plt.legend(fontsize=14)
plt.subplot(3, 1, 3)
sns.kdeplot(np.array(posteriors_DSAS["S0"]), color='black', fill=True, alpha=1)
x = np.linspace(-12.01, -2.99, 1000)
plt.plot(x, uniform.pdf(x, loc=-12, scale=(-3--12)), color="slategrey")
plt.fill_between(x=[-12, -2.99], y1=uniform.pdf(x[10], loc=-12, scale=(-3--12)), y2=0, color="slategrey", alpha=0.6)
plt.xlabel('log(S0)', fontsize=14)
plt.yticks([0, 1], labels=['',''])
plt.ylabel('')
plt.tight_layout()
plt.show()
plt.savefig('KLD_histograms.png', dpi=300)
plt.close()
