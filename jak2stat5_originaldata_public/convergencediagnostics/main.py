import numpy as np
from convergence_diagnostics import convergencediagonistics_logscale
from split_r_convergence_metric import split_r_convergence_metric
from split_r_convergence_metric import qqplot

#master variables
proposal_notation = 'proposal 0.1'
data_notation = 'baseline'

#convergence diagnostics and metrics
chain1 = np.load('generateinputs/chain1_{}_{}_estimatedvariables.npy'.format(data_notation, proposal_notation))
chain2 = np.load('generateinputs/chain2_{}_{}_estimatedvariables.npy'.format(data_notation, proposal_notation))
chain3 = np.load('generateinputs/chain3_{}_{}_estimatedvariables.npy'.format(data_notation, proposal_notation))
chain4 = np.load('generateinputs/chain4_{}_{}_estimatedvariables.npy'.format(data_notation, proposal_notation))
chain5 = np.load('generateinputs/chain5_{}_{}_estimatedvariables.npy'.format(data_notation, proposal_notation))

#generate supplementary figure S7A
convergencediagonistics_logscale(chain1, chain2, chain3, chain4, chain5)

#generate supplementary figure6A
split_r_convergence_metric(chain1, chain2, chain3, chain4, chain5)

#generate Q-Q plot to check Gelman Rubin applicability, figure in Section 1.3 of Supplement
qqplot(np.load('generateinputs/all_chains_{}_{}_estimatedvariables.npy'.format(data_notation, proposal_notation)))

