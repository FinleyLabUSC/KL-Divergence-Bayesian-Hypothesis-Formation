from STAT_JAK_BCL_forward_uncertainty_propogation import STAT_JAK_BCL_forward_uncertainty_propogation
import numpy as np

#generate supplementary figure S9
simulations = np.load('generateinputs/fine_grain_JAK_STAT_BCL_simulation.npy', allow_pickle=True)
STAT_JAK_BCL_forward_uncertainty_propogation(simulations, 'extended data', 'slategrey', 0.05, 0.95 )