from shuffle import isolate_estimated_variables
import numpy as np
from scipy.integrate import solve_ivp
from odemodels import simulatejakstat
import calculate


# master variables
proposal_notation = 'proposal 0.1'
data_notation = 'baseline and timepoint'

chain1 = np.load('generateinputs/chain1_{}_{}_estimatedvariables.npy'.format(data_notation, proposal_notation))
chain2 = np.load('generateinputs/chain2_{}_{}_estimatedvariables.npy'.format(data_notation, proposal_notation))
chain3 = np.load('generateinputs/chain3_{}_{}_estimatedvariables.npy'.format(data_notation, proposal_notation))
totalchain = np.concatenate((chain1, chain2, chain3), axis=0)
np.save('generateinputs/all_chains_{}_{}_estimatedvariables.npy'.format(data_notation, proposal_notation), totalchain)

#subsample estimated variables to ~10000 samples
totalchain = np.load('generateinputs/all_chains_{}_{}_estimatedvariables.npy'.format(data_notation, proposal_notation))
thinned = totalchain[0:-1:int((len(totalchain)/10000)), :]
np.save('generateinputs/all_chains_{}_{}_thinnedestimatedvariables.npy'.format(data_notation, proposal_notation), thinned)

# generate fine grain simulation with subsampled samples
parameters_subsampled = np.load('generateinputs/all_chains_{}_{}_thinnedparameters.npy'.format(data_notation, proposal_notation))
species_subsampled = np.load('generateinputs/all_chains_{}_{}_thinnedspecies.npy'.format(data_notation, proposal_notation))

tspan_fine = np.load('generateinputs/tspan_fine_all_sec.npy')
species_dictionary = np.load('generateinputs/dictionary_all_species.npy')
pstata = []
pstatb = []
transloc_preda = []
transloc_predb = []
bcl = []
pjak = []
socs = []

len(species_subsampled[:, 1])

for i in np.arange(0, len(species_subsampled[:, 1])):
    print(i)
    species = species_subsampled[i, :]
    parameters = parameters_subsampled[i, :]
    sol = solve_ivp(simulatejakstat, [0, np.max(tspan_fine)], species.flatten(), method='BDF',
                    args=(parameters.flatten(),), t_eval=tspan_fine, rtol=1e-9, atol=1e-12)
    simulation = np.transpose(sol.y)
    SOCSindex1 = np.where(species_dictionary == 'SOCS1')[0][0]
    SOCSindex2 = np.where(species_dictionary == 'SOCS1PRLRJ2a')[0][0]
    SOCSindex3 = np.where(species_dictionary == 'SOCS1PRLRJ2aSHP2')[0][0]
    total_SOCS = simulation[:, SOCSindex1] + simulation[:, SOCSindex2] + simulation[:, SOCSindex3]
    normalized_SOCS = np.divide(total_SOCS, np.max(total_SOCS))
    socs.append(normalized_SOCS)
    pStatA_norm, pStatB_norm, transloc_predA, transloc_predB, BcL_foldchange, pJAK2_norm = calculate.experimental_predictions_fine_tspan(simulation)
    pstata.append(pStatA_norm)
    pstatb.append(pStatB_norm)
    transloc_preda.append(transloc_predA)
    transloc_predb.append(transloc_predB)
    pjak.append(pJAK2_norm)
    bcl.append(BcL_foldchange)
fine_grain_simulation = np.array((pstata, pstatb, transloc_preda, transloc_predb, pjak, bcl), dtype=object)
socs_fine_grain_simulation = np.array(socs)
np.save('generateinputs/fine_grain_JAK_STAT_BCL_simulation.npy', fine_grain_simulation)
np.save('generateinputs/fine_grain_SOCS_simulation.npy', socs_fine_grain_simulation)

#generate k24 socs sorted samples
k24_lower_sample_cutoff = 0.3
k24_upper_sample_cutoff = 0.6
socsimulation = np.load('generateinputs/fine_grain_SOCS_simulation.npy')
posteriorsamples_thinned = np.load('generateinputs/all_chains_{}_{}_thinnedestimatedvariables.npy'.format(data_notation, proposal_notation))
variable_dictionary = np.load('generateinputs/dictionary_estimatedvariables.npy')
k24_index = np.where(variable_dictionary == 'k24')[0][0]
k24_values = posteriorsamples_thinned[:, k24_index]
matrix = np.transpose(np.concatenate(([k24_values], np.transpose(socsimulation)), axis=0))

posteriorsamples_thinned_sorted = matrix[k24_values.argsort()]
lower_index = int(k24_lower_sample_cutoff * len(posteriorsamples_thinned_sorted))
upper_index = int(k24_upper_sample_cutoff * len(posteriorsamples_thinned_sorted))
max_k24 = posteriorsamples_thinned_sorted[-upper_index:-1:2, 1:4323]
min_k24 = posteriorsamples_thinned_sorted[0:lower_index, 1:4323]

np.save('generateinputs/fine_grain_SOCS_simulation_sorted_by_upper_{}%_k24'.format(k24_upper_sample_cutoff*100), max_k24)
np.save('generateinputs/fine_grain_SOCS_simulation_sorted_by_lower_{}%_k24'.format(k24_lower_sample_cutoff*100), min_k24)

#generate k24 experimental sorted samples
k24_lower_sample_cutoff = 0.3
k24_upper_sample_cutoff = 0.6

experimental_proteins = np.load('generateinputs/fine_grain_JAK_STAT_BCL_simulation.npy', allow_pickle=True)
posteriorsamples_thinned = np.load('generateinputs/all_chains_{}_{}_thinnedestimatedvariables.npy'.format(data_notation, proposal_notation))
variable_dictionary = np.load('generateinputs/dictionary_estimatedvariables.npy')
k24_index = np.where(variable_dictionary == 'k24')[0][0]
k24_values = posteriorsamples_thinned[:, k24_index]
sorted_samples = experimental_proteins[:, k24_values.argsort(), :]

lower_index = int(k24_lower_sample_cutoff * len(sorted_samples[1, :, 1]))
upper_index = int(k24_upper_sample_cutoff * len(sorted_samples[1, :, 1]))
max_k24 = sorted_samples[:, -upper_index:-1:2, :]
min_k24 = sorted_samples[:, 0:lower_index, :]

np.save('generateinputs/fine_grain_JAK_STAT_BCL_simulation_sorted_by_upper_{}%_k24'.format(k24_upper_sample_cutoff*100), max_k24)
np.save('generateinputs/fine_grain_JAK_STAT_BCL_simulation_sorted_by_lower_{}%_k24'.format(k24_lower_sample_cutoff*100), min_k24)

#store k24 samples
k24_lower_sample_cutoff = 0.3
k24_upper_sample_cutoff = 0.6

posteriorsamples_thinned = np.load('generateinputs/all_chains_{}_{}_thinnedestimatedvariables.npy'.format(data_notation, proposal_notation))
variable_dictionary = np.load('generateinputs/dictionary_estimatedvariables.npy')
k24_index = np.where(variable_dictionary == 'k24')[0][0]
k24_values = posteriorsamples_thinned[:, k24_index]
sorted_samples = k24_values[k24_values.argsort()]
lower_index = int(k24_lower_sample_cutoff * len(sorted_samples))
upper_index = int(k24_upper_sample_cutoff * len(sorted_samples))
max_k24 = sorted_samples[-upper_index:-1:2]
min_k24 = sorted_samples[0:lower_index]

np.save('generateinputs/k24_samples_sorted_by_upper_{}%_k24'.format(k24_upper_sample_cutoff*100), max_k24)
np.save('generateinputs/k24_samples_sorted_by_lower_{}%_k24'.format(k24_lower_sample_cutoff*100), min_k24)
np.save('generateinputs/k24_samples_all', k24_values)