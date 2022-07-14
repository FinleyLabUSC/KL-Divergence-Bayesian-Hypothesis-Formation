import numpy as np
import matplotlib.pyplot as plt

def STAT_JAK_BCL_forward_uncertainty_propogation(simulations, case, color, lower_percentile,upper_percentile):
        tspan_fine = np.load('modelverification/tspan_fine_all_sec.npy')
        tspan_course = np.load('modelverification/tspan_course_all_sec.npy')
        pstata = simulations[0, :, :]
        pstatb = simulations[1, :, :]
        transloca = simulations[2, :, :]
        translocb = simulations[3, :, :]
        pjak = simulations[4, :, :]
        bcl = simulations[5, :, :]
        exp_data = np.load('modelverification/data_experimental.npy').flatten()
        exp_mu_bcl = exp_data[26:33]
        exp_mu_pjak = exp_data[33:39]
        exp_mu_pstata = exp_data[0:6]
        exp_mu_pstatb = exp_data[6:12]
        exp_mu_transloca = exp_data[12:17]
        exp_mu_translocb = exp_data[17:26]

        exp_tspan_pjak_pstat = np.around(np.divide(tspan_course[[2, 4, 5, 7, 9, 11]], 3600), 2)
        exp_tspan_transloca = np.around(np.divide(tspan_course[[4, 5, 6, 8, 11]], 3600), 2)
        exp_tspan_translocb = np.around(np.divide(tspan_course[[1, 3, 4, 5, 7, 8, 9, 10, 11]], 3600), 2)
        exp_tspan_bcl = np.around(np.divide(tspan_course[[7, 9, 11, 12, 13, 14, 15]], 3600), 2)

        sim = [pstata, transloca, pstatb, translocb, pjak, bcl]
        exp = [exp_mu_pstata, exp_mu_transloca, exp_mu_pstatb, exp_mu_translocb, exp_mu_pjak, exp_mu_bcl]
        names = ['[pSTATA]', 'N:C [STATA]', '[pSTATB]', 'N:C [STATB]', '[pJAK]', '[BCL-xL]']
        exp_tspan = [exp_tspan_pjak_pstat, exp_tspan_transloca, exp_tspan_pjak_pstat, exp_tspan_translocb,
                     exp_tspan_pjak_pstat, exp_tspan_bcl]

        plt.figure(figsize=(10, 7))
        for s, e, t, n, j in zip(sim, exp, exp_tspan, names, np.arange(1, 7)):
            plt.subplot(3, 2, j)
            plt.subplots_adjust(bottom=0.1, top=0.99, left=0.05, right=0.95)
            quantiles_low = np.array(np.quantile(s, lower_percentile, axis=0), dtype=float)
            quantiles_high = np.array(np.quantile(s, upper_percentile, axis=0), dtype=float)
            median = np.array(np.median(s, axis=0), dtype=float)
            plt.fill_between(tspan_fine, quantiles_low, quantiles_high, alpha=0.3, color=color,
                             label='{}th to {}th percentile'.format(lower_percentile * 100, upper_percentile * 100))
            plt.plot(tspan_fine, median, color=color, linestyle='dotted', label='median prediction', linewidth=3)
            plt.plot(t * 3600, e.flatten(), 'o', color='black', label='experimental data mean')
            if n == '[pJAK]':
                plt.xlabel('Time, hours', fontsize=14)
            if n == '[BCL-xL]':
                plt.xlabel('Time, hours', fontsize=14)
            # if n == '[pSTATB]':
            # plt.legend(loc='upper right', fontsize=8)
            plt.ylim([0, 3.5])
            plt.xlim([0, t.max() * 3600 + 1000])
            plt.xticks(fontsize=8, ticks=t * 3600, labels=t, rotation=90)
            if n == '[pSTATA]' or n == '[pSTATB]':
                plt.xticks(fontsize=8, ticks=t * 3600, labels=['','','','','',''], rotation=90)
            plt.yticks(fontsize=8)
            plt.ylabel('relative {}'.format(n), fontsize=14)
        plt.savefig('modelverification/JAK_STAT_BCL_timecourses_{}_{}th_to_{}th_percentile_edited.png'.format(case,
                                                                                                       lower_percentile * 100,
                                                                                                       upper_percentile * 100), dpi=300)
        plt.show()