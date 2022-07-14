from metropolishastings import metropolishastingsalgo
from metropolishastings import initialize_at_SOCS_modes
from metropolishastings import initializeMH
import numpy as np
from likelihoods import log_likelihood_baseline_normal_inverse_gamma_variance
from datetime import datetime
from shuffle import isolate_estimated_variables

#MCMC sampler - samples analyzed in manuscript are saved for analysis, so there is no need to rerun sampler.
# However, code is made available for cross-checking and may be used for future projects.

#master variables
likelihoodfunction = log_likelihood_baseline_normal_inverse_gamma_variance
prior = np.sqrt(2)
proposal = 0.1
proposal_notation = 'proposal 0.1'
data_notation = 'baseline'

#generate metropolois hasting initializations from random start
initializeMH(np.load('posteriorapproximation/initialization_matrix_prior_guess_.npy', allow_pickle=True),
             prior, likelihoodfunction)

# generate metropolis hastings initializations from k24 modes
max_k24 = np.load('posteriorapproximation/initialization_maxk24.npy')
min_k24 = np.load('posteriorapproximation/initialization_mink24.npy')
initialize_at_SOCS_modes(max_k24, min_k24, likelihoodfunction, 0)

# metropolis hastings algorithm
#setting 0
initialize_chain1 = np.load('posteriorapproximation/initialization_matrix_prior_guess_.npy', allow_pickle=True)
initialize_chain2 = np.load('posteriorapproximation/initialization_matrix_minlikelihood.npy', allow_pickle=True)
initialize_chain3 = np.load('posteriorapproximation/initialization_matrix_maxlikelihood.npy', allow_pickle=True)
initialize_chain4 = np.load('posteriorapproximation/initialization_matrix_min k24.npy', allow_pickle=True)
initialize_chain5 = np.load('posteriorapproximation/initialization_matrix_max k24.npy', allow_pickle=True)

#setting 1
initialize_chain1 = np.array((np.load('posteriorapproximation/chain1_{}_{}_parameters.npy'.format(data_notation, proposal_notation)),
                             np.load('posteriorapproximation/chain1_{}_{}_species.npy'.format(data_notation, proposal_notation)),
                             np.load('posteriorapproximation/chain1_{}_{}_likelihood.npy'.format(data_notation, proposal_notation)),
                             np.load('posteriorapproximation/chain1_{}_{}_predictions.npy'.format(data_notation, proposal_notation)),
                             np.load('posteriorapproximation/chain1_{}_{}_acceptancerate.npy'.format(data_notation, proposal_notation))), dtype=object)
initialize_chain2 = np.array((np.load('posteriorapproximation/chain2_{}_{}_parameters.npy'.format(data_notation, proposal_notation)),
                             np.load('posteriorapproximation/chain2_{}_{}_species.npy'.format(data_notation, proposal_notation)),
                             np.load('posteriorapproximation/chain2_{}_{}_likelihood.npy'.format(data_notation, proposal_notation)),
                             np.load('posteriorapproximation/chain2_{}_{}_predictions.npy'.format(data_notation, proposal_notation)),
                             np.load('posteriorapproximation/chain2_{}_{}_acceptancerate.npy'.format(data_notation, proposal_notation))), dtype=object)
initialize_chain3 = np.array((np.load('posteriorapproximation/chain3_{}_{}_parameters.npy'.format(data_notation, proposal_notation)),
                             np.load('posteriorapproximation/chain3_{}_{}_species.npy'.format(data_notation, proposal_notation)),
                             np.load('posteriorapproximation/chain3_{}_{}_likelihood.npy'.format(data_notation, proposal_notation)),
                             np.load('posteriorapproximation/chain3_{}_{}_predictions.npy'.format(data_notation, proposal_notation)),
                             np.load('posteriorapproximation/chain3_{}_{}_acceptancerate.npy'.format(data_notation, proposal_notation))), dtype=object)
initialize_chain4 = np.array((np.load('posteriorapproximation/chain4_{}_{}_parameters.npy'.format(data_notation, proposal_notation)),
                             np.load('posteriorapproximation/chain4_{}_{}_species.npy'.format(data_notation, proposal_notation)),
                             np.load('posteriorapproximation/chain4_{}_{}_likelihood.npy'.format(data_notation, proposal_notation)),
                             np.load('posteriorapproximation/chain4_{}_{}_predictions.npy'.format(data_notation, proposal_notation)),
                             np.load('posteriorapproximation/chain4_{}_{}_acceptancerate.npy'.format(data_notation, proposal_notation))), dtype=object)
initialize_chain5 = np.array((np.load('posteriorapproximation/chain5_{}_{}_parameters.npy'.format(data_notation, proposal_notation)),
                             np.load('posteriorapproximation/chain5_{}_{}_species.npy'.format(data_notation, proposal_notation)),
                             np.load('posteriorapproximation/chain5_{}_{}_likelihood.npy'.format(data_notation, proposal_notation)),
                             np.load('posteriorapproximation/chain5_{}_{}_predictions.npy'.format(data_notation, proposal_notation)),
                             np.load('posteriorapproximation/chain5_{}_{}_acceptancerate.npy'.format(data_notation, proposal_notation))), dtype=object)

initialize = [initialize_chain1, initialize_chain2, initialize_chain3, initialize_chain4, initialize_chain5]
path1 = ['posteriorapproximation/chain1_{}_{}_parameters.npy'.format(data_notation, proposal_notation), 'posteriorapproximation/chain1_{}_{}_species.npy'.format(data_notation, proposal_notation),
         'posteriorapproximation/chain1_{}_{}_likelihood.npy'.format(data_notation, proposal_notation), 'posteriorapproximation/chain1_{}_{}_predictions.npy'.format(data_notation, proposal_notation),
         'posteriorapproximation/chain1_{}_{}_simulations.npy'.format(data_notation, proposal_notation), 'posteriorapproximation/chain1_{}_{}_acceptancerate.npy'.format(data_notation, proposal_notation)]
path2 = ['posteriorapproximation/chain2_{}_{}_parameters.npy'.format(data_notation, proposal_notation), 'posteriorapproximation/chain2_{}_{}_species.npy'.format(data_notation, proposal_notation),
         'posteriorapproximation/chain2_{}_{}_likelihood.npy'.format(data_notation, proposal_notation), 'posteriorapproximation/chain2_{}_{}_predictions.npy'.format(data_notation, proposal_notation),
         'posteriorapproximation/chain2_{}_{}_simulations.npy'.format(data_notation, proposal_notation), 'posteriorapproximation/chain2_{}_{}_acceptancerate.npy'.format(data_notation, proposal_notation)]
path3 = ['posteriorapproximation/chain3_{}_{}_parameters.npy'.format(data_notation, proposal_notation), 'posteriorapproximation/chain3_{}_{}_species.npy'.format(data_notation, proposal_notation),
         'posteriorapproximation/chain3_{}_{}_likelihood.npy'.format(data_notation, proposal_notation), 'posteriorapproximation/chain3_{}_{}_predictions.npy'.format(data_notation, proposal_notation),
         'posteriorapproximation/chain3_{}_{}_simulations.npy'.format(data_notation, proposal_notation), 'posteriorapproximation/chain3_{}_{}_acceptancerate.npy'.format(data_notation, proposal_notation)]
path4 = ['posteriorapproximation/chain4_{}_{}_parameters.npy'.format(data_notation, proposal_notation), 'posteriorapproximation/chain4_{}_{}_species.npy'.format(data_notation, proposal_notation),
         'posteriorapproximation/chain4_{}_{}_likelihood.npy'.format(data_notation, proposal_notation), 'posteriorapproximation/chain4_{}_{}_predictions.npy'.format(data_notation, proposal_notation),
         'posteriorapproximation/chain4_{}_{}_simulations.npy'.format(data_notation, proposal_notation), 'posteriorapproximation/chain4_{}_{}_acceptancerate.npy'.format(data_notation, proposal_notation)]
path5 = ['posteriorapproximation/chain5_{}_{}_parameters.npy'.format(data_notation, proposal_notation), 'posteriorapproximation/chain5_{}_{}_species.npy'.format(data_notation, proposal_notation),
         'posteriorapproximation/chain5_{}_{}_likelihood.npy'.format(data_notation, proposal_notation), 'posteriorapproximation/chain5_{}_{}_predictions.npy'.format(data_notation, proposal_notation),
         'posteriorapproximation/chain5_{}_{}_simulations.npy'.format(data_notation, proposal_notation), 'posteriorapproximation/chain5_{}_{}_acceptancerate.npy'.format(data_notation, proposal_notation)]

paths = [path1, path2, path3, path4, path5]

# MCMC
init_time = datetime.now()
for i in np.arange(0, 5):
    initialize_chain = initialize[i]
    pathi = paths[i]
    prior_std = prior
    proposal_std = proposal
    likelihood_sc, biochemical_species, biochemical_parameters, acceptance_rate, exp_predictions = metropolishastingsalgo(initialize_chain, 100000, prior_std, proposal_std, likelihoodfunction, 1)
    np.save(pathi[0], np.array(biochemical_parameters))
    np.save(pathi[1], np.array(biochemical_species))
    np.save(pathi[2], np.array(likelihood_sc))
    np.save(pathi[3], np.array(exp_predictions))
    # np.save(pathi[4], np.array(simulations))
    np.save(pathi[5], np.array(acceptance_rate))
fin_time = datetime.now()
total_time = fin_time-init_time
print('total time:', total_time)


#now that samples have been generated, we process them for compatability with future analyses:

# name variables
proposal_notation = 'proposal 0.1'
data_notation = 'baseline'

# isolate estimated variables and save in generateinputs folder
inputs = [['posteriorapproximation/chain1_{}_{}_parameters.npy'.format(data_notation, proposal_notation), 'posteriorapproximation/chain1_{}_{}_species.npy'.format(data_notation, proposal_notation),
          'generateinputs/chain1_{}_{}_estimatedvariables.npy'.format(data_notation, proposal_notation)],
          ['posteriorapproximation/chain2_{}_{}_parameters.npy'.format(data_notation, proposal_notation), 'posteriorapproximation/chain2_{}_{}_species.npy'.format(data_notation, proposal_notation),
          'generateinputs/chain2_{}_{}_estimatedvariables.npy'.format(data_notation, proposal_notation)],
          ['posteriorapproximation/chain3_{}_{}_parameters.npy'.format(data_notation, proposal_notation), 'posteriorapproximation/chain3_{}_{}_species.npy'.format(data_notation, proposal_notation),
          'generateinputs/chain3_{}_{}_estimatedvariables.npy'.format(data_notation, proposal_notation)],
          ['posteriorapproximation/chain4_{}_{}_parameters.npy'.format(data_notation, proposal_notation), 'posteriorapproximation/chain4_{}_{}_species.npy'.format(data_notation, proposal_notation),
          'generateinputs/chain4_{}_{}_estimatedvariables.npy'.format(data_notation, proposal_notation)],
          ['posteriorapproximation/chain5_{}_{}_parameters.npy'.format(data_notation, proposal_notation), 'posteriorapproximation/chain5_{}_{}_species.npy'.format(data_notation, proposal_notation),
          'generateinputs/chain5_{}_{}_estimatedvariables.npy'.format(data_notation, proposal_notation)]]

for i in np.arange(0, 5):
    inputi = inputs[i]
    isolate_estimated_variables(np.load(inputi[0]), np.load(inputi[1]), inputi[2],)

# thin parameters and species matrices and save in generateinputs folder

#subsample parameters and species to ~10000 samples
parameters = np.concatenate((np.load('posteriorapproximation/chain1_{}_{}_parameters.npy'.format(data_notation, proposal_notation))[1000:-1, :],
                             np.load('posteriorapproximation/chain2_{}_{}_parameters.npy'.format(data_notation, proposal_notation))[1000:-1, :],
                              np.load('posteriorapproximation/chain3_{}_{}_parameters.npy'.format(data_notation, proposal_notation))[1000:-1, :],
                               np.load('posteriorapproximation/chain4_{}_{}_parameters.npy'.format(data_notation, proposal_notation))[1000:-1, :],
                                np.load('posteriorapproximation/chain5_{}_{}_parameters.npy'.format(data_notation, proposal_notation))[1000:-1, :]))

species = np.concatenate((np.load('posteriorapproximation/chain1_{}_{}_species.npy'.format(data_notation, proposal_notation))[1000:-1, :],
                             np.load('posteriorapproximation/chain2_{}_{}_species.npy'.format(data_notation, proposal_notation))[1000:-1, :],
                              np.load('posteriorapproximation/chain3_{}_{}_species.npy'.format(data_notation, proposal_notation))[1000:-1, :],
                               np.load('posteriorapproximation/chain4_{}_{}_species.npy'.format(data_notation, proposal_notation))[1000:-1, :],
                                np.load('posteriorapproximation/chain5_{}_{}_species.npy'.format(data_notation, proposal_notation))[1000:-1, :]))

parameters_thinned = parameters[0:-1:int((len(parameters)/10000)), :]
species_thinned = species[0:-1:int((len(species)/10000)), :]
np.save('generateinputs/all_chains_{}_{}_thinnedparameters.npy'.format(data_notation, proposal_notation), parameters_thinned)
np.save('generateinputs/all_chains_{}_{}_thinnedspecies.npy'.format(data_notation, proposal_notation), species_thinned)

