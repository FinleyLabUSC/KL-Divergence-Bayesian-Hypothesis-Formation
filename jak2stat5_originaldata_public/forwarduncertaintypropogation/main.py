from STAT_JAK_BCL_k24_forward_uncertainty_propogation import STAT_JAK_BCL_k24_forward_uncertainty_propogation
import numpy as np
from k24_samples import visualize_k24_samples
from SOCS_k24_forward_uncertainty_propogation import SOCS_k24_forward_uncertainty_propogation

#generate supplementary figure S8
upperk24simulations = np.load('generateinputs/fine_grain_JAK_STAT_BCL_simulation_sorted_by_upper_60.0%_k24.npy', allow_pickle=True)
lowerk24simulations = np.load('generateinputs/fine_grain_JAK_STAT_BCL_simulation_sorted_by_lower_30.0%_k24.npy', allow_pickle=True)
simulations = [upperk24simulations, lowerk24simulations]
color = ['orange', 'darkturquoise']
k24 = ['upper', 'lower']

for i in np.arange(0, 2):
    STAT_JAK_BCL_k24_forward_uncertainty_propogation(simulations[i], 'forwarduncertaintypropogation', color[i],
                                                     0.25, 0.75, k24[i])

#generate figure 4A
allk24 = np.load('generateinputs/k24_samples_all.npy')
lowerk24 = np.load('generateinputs/k24_samples_sorted_by_lower_30.0%_k24.npy')
upperk24 = np.load('generateinputs/k24_samples_sorted_by_upper_60.0%_k24.npy')
visualize_k24_samples(allk24, upperk24, lowerk24,'forwarduncertaintypropogation')

#generate figure 4B
upperk24SOCsimulations = np.load('generateinputs/fine_grain_SOCS_simulation_sorted_by_upper_60.0%_k24.npy')
lowerk24SOCsimulations = np.load('generateinputs/fine_grain_SOCS_simulation_sorted_by_lower_30.0%_k24.npy')
SOCS_k24_forward_uncertainty_propogation(upperk24SOCsimulations, lowerk24SOCsimulations,
                                         'forwarduncertaintypropogation', 0.25, 0.75, 18 )