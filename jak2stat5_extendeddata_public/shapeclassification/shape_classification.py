import numpy as np
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
import squarify

def shape_classification(socs_simulations):
    no_peak_and_sustained = []
    no_peak = []
    sustained = []
    early_peak_and_transient = []
    indicator = []
    socs_simulations_4hrs = socs_simulations[:, 0:int(4*3600/20)]
    for i in np.arange(0, len(socs_simulations[:, 1])):
        peaks_index = find_peaks(socs_simulations_4hrs[i, :])[0]
        if peaks_index.size == 0 and socs_simulations[i, int(24*3600/20)] > 0.5:
            indicator.append(0)
            no_peak_and_sustained.append(socs_simulations[i, :])
        elif peaks_index.size == 0:
            indicator.append(0)
            no_peak.append(socs_simulations[i, :])
        elif socs_simulations[i, int(24*3600/20)] > 0.5:
            indicator.append(0)
            sustained.append(socs_simulations[i, :])
        elif peaks_index.size != 0 and socs_simulations[i, int(24*3600/20)] <= 0.5:
            early_peak_and_transient.append(socs_simulations[i, :])
            indicator.append(1)
    check = len(no_peak_and_sustained)+len(no_peak)+len(sustained)+len(early_peak_and_transient)-len(socs_simulations[:, 1])
    if check == 0:
        print('All simulations accounted for :)')
    if check != 0:
        print('Shape classification failed :/')

    simulations = [np.array(early_peak_and_transient).transpose(), np.array(no_peak_and_sustained).transpose(),
                   np.array(no_peak).transpose(), np.array(sustained).transpose()]
    colors = ['darkorchid', 'black', 'black', 'black']
    shape = ['early peak, transient', 'no peak, sustained', 'no peak, transient', 'early peak, sustained']

    for i in np.arange(0, 4):
        sim = simulations[i]
        c = colors[i]
        s = shape[i]
        tspan_course = np.load('shapeclassification/tspan_course_all_sec.npy')
        plt.plot(sim[:, 0:-1:10], alpha=0.1, color=c)
        plt.subplots_adjust(top=0.99, right=0.99)
        plt.xticks(ticks=tspan_course[[5, 7, 8, 9, 10, 11, 12, 13, 14, 15]] / 20, labels=np.divide(tspan_course[[5, 7, 8, 9, 10, 11, 12, 13, 14, 15]], 3600).astype(int))
        plt.xlabel('time, hours')
        plt.ylabel('relative [SOCS]')
        plt.savefig('shapeclassification/line_plot_{}.png'.format(s), dpi=300)
        plt.close()

    size = [len(early_peak_and_transient), len(no_peak_and_sustained), len(no_peak), len(sustained)]
    percentage = np.around(np.divide(size, len(socs_simulations[:, 1]))*100, decimals=0)
    colors2 = ['darkorchid', 'silver', 'whitesmoke', 'dimgrey']
    shape= ['early peak, transient:\n {}%'.format(int(percentage[0])), 'no peak, sustained: {}%'.format(int(percentage[1])), 'no peak, transient:\n {}%'.format(int(percentage[2])), 'early peak, sustained:\n {}%'.format(int(percentage[3]))]
    plt.rc('font', size=14)
    squarify.plot(sizes=size, label=shape, color=colors2)
    plt.axis('off')
    plt.savefig('shapeclassification/treemap_plot.png', dpi=300)

    np.save('shapeclassification/literature_shape_indicator.npy', indicator)

def adjustable_shape_classification(socs_simulations, short_term_criterion, long_term_criterion):
    no_peak_and_sustained = []
    no_peak = []
    sustained = []
    early_peak_and_transient = []
    indicator = []
    socs_simulations_4hrs = socs_simulations[:, 0:int(short_term_criterion*3600/20)]
    for i in np.arange(0, len(socs_simulations[:, 1])):
        peaks_index = find_peaks(socs_simulations_4hrs[i, :])[0]
        if peaks_index.size == 0 and socs_simulations[i, int(24*3600/20)] > long_term_criterion:
            indicator.append(0)
            no_peak_and_sustained.append(socs_simulations[i, :])
        elif peaks_index.size == 0:
            indicator.append(0)
            no_peak.append(socs_simulations[i, :])
        elif socs_simulations[i, int(24*3600/20)] > long_term_criterion:
            indicator.append(0)
            sustained.append(socs_simulations[i, :])
        elif peaks_index.size != 0 and socs_simulations[i, int(24*3600/20)] <= long_term_criterion:
            early_peak_and_transient.append(socs_simulations[i, :])
            indicator.append(1)
    check = len(no_peak_and_sustained)+len(no_peak)+len(sustained)+len(early_peak_and_transient)-len(socs_simulations[:, 1])
    if check == 0:
        print('All simulations accounted for :)')
    if check != 0:
        print('Shape classification failed :/')

    size = [len(early_peak_and_transient), len(no_peak_and_sustained), len(no_peak), len(sustained)]
    percentage = np.around(np.divide(size, len(socs_simulations[:, 1]))*100, decimals=0)
    colors2 = ['darkorchid', 'silver', 'whitesmoke', 'dimgrey']
    shape= ['{}%'.format(int(percentage[0])), '{}%'.format(int(percentage[1])), '{}%'.format(int(percentage[2])), '{}%'.format(int(percentage[3]))]
    plt.rc('font', size=14)
    squarify.plot(sizes=size, label=shape, color=colors2)
    plt.axis('off')
    plt.savefig('shapeclassification/treemap_plot_{}hrs_{}threshold.png'.format(short_term_criterion, long_term_criterion), dpi=300)





