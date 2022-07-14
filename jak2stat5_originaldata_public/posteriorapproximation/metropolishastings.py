import numpy as np
import shuffle
import calculate
from prior_and_proposal import generate_lognormal
from prior_and_proposal import evaluate_lognormal_function
from scipy.integrate import solve_ivp
from odemodels import simulatejakstat
import math
import matplotlib.pyplot as plt

def initializeMH(initialization, prior_std, likelihood):
    #initialize arrays
    parameters = initialization[0]
    species = initialization[1]
    likelihood_sc = initialization[2]
    exp_pred = initialization[3]
    sim = initialization[4]

    # initialize parameter, species, likelihood, simulation, acceptance rate lists to save
    biochemical_parameters = [parameters, parameters]
    biochemical_species = [species, species]
    likelihood_score = [likelihood_sc, likelihood_sc]
    exp_predictions = [exp_pred, exp_pred]
    simulations = [sim, sim]

    # initialize initial guess for parameter&species matrix
    fixed_params = np.load('posteriorapproximation/prior_parameters.npy')
    fixed_species = np.load('posteriorapproximation/prior_species.npy')

    #initialize tspan for ODE simulation
    tspan_course = np.load('posteriorapproximation/tspan_course_all_sec.npy')
    tspan_fine = np.load('posteriorapproximation/tspan_fine_all_sec.npy')
    indices = []
    for i in tspan_course:
        index = np.argwhere(tspan_fine == i)
        indices.append(index)
    indices = np.array(indices).flatten()

    for i in np.arange(0, 10000):
        print(i)
        # generate proposed parameterization
        estimated_parameters = shuffle.index_estimated_parameters(fixed_params)
        estimated_species = shuffle.index_estimated_species(fixed_species)
        sampled_estimated_parameters = generate_lognormal(estimated_parameters, prior_std)
        sampled_estimated_species = generate_lognormal(estimated_species, prior_std)
        proposed_parameters = shuffle.assign_estimated_parameters(sampled_estimated_parameters, fixed_params)
        proposed_species = shuffle.assign_estimated_species(sampled_estimated_species, fixed_species)
        proposed_parameters, proposed_species = calculate.dependent_values(proposed_species, proposed_parameters)

        # simulate with proposed parameters
        sol = solve_ivp(simulatejakstat, [0, np.max(tspan_fine)], proposed_species.flatten(), method='BDF',
                        args=(proposed_parameters.flatten(),), t_eval=tspan_fine, rtol=1e-9, atol=1e-12)
        proposed_simulation_fine = np.transpose(sol.y)
        proposed_simulation = proposed_simulation_fine[indices, :]

        # calculate experimental predictions
        pStatA_norm_prop, pStatB_norm_prop, transloc_predA_prop, transloc_predB_prop, BcL_foldchange_prop, pJAK2_norm_prop\
            = calculate.experimental_predictions(proposed_simulation)
        proposed_predictions = np.concatenate((BcL_foldchange_prop, pJAK2_norm_prop,pStatA_norm_prop,
                                               pStatB_norm_prop, transloc_predA_prop, transloc_predB_prop), axis=None)
        # log likelihood
        log_likelihood_proposed = likelihood(proposed_predictions)
        print('log likelihood', log_likelihood_proposed)
        # append
        biochemical_parameters.append(proposed_parameters)
        biochemical_species.append(proposed_species)
        likelihood_score.append(log_likelihood_proposed)
        exp_predictions.append(proposed_predictions)
        simulations.append(proposed_simulation_fine)
    indexmax = np.argmax(likelihood_score)
    indexmin = np.argmin(likelihood_score)

    # format and save initialization matrix
    initialization_max = np.array([biochemical_parameters[indexmax], biochemical_species[indexmax], likelihood_score[indexmax], exp_predictions[indexmax], simulations[indexmax]], dtype=object)
    initialization_min = np.array([biochemical_parameters[indexmin], biochemical_species[indexmin], likelihood_score[indexmin], exp_predictions[indexmin], simulations[indexmin]], dtype=object)
    np.save('posteriorapproximation/initialization_matrix_maxlikelihood.npy', initialization_max)
    np.save('posteriorapproximation/initialization_matrix_minlikelihood.npy', initialization_min)

def initialize_at_SOCS_modes(max_k24, min_k24, likelihood, setting):
    # function plots SOCS timecourse and generates initialization matrix
    # note - plot does not change between simulations, but initialization of likelihood does
    # setting determines whether plot is generated

    species_dictionary = np.load('posteriorapproximation/dictionary_all_species.npy')

    # initialize tspan for ODE simulation
    tspan_course = np.load('posteriorapproximation/tspan_course_all_sec.npy')
    tspan_fine = np.load('posteriorapproximation/tspan_fine_all_sec.npy')
    indices = []
    for i in tspan_course:
        index = np.argwhere(tspan_fine == i)
        indices.append(index)
    indices = np.array(indices).flatten()

    prior_parameters = np.load('posteriorapproximation/prior_parameters.npy')
    prior_species = np.load('posteriorapproximation/prior_species.npy')

    color = ['orange', 'darkturquoise']
    variables = [max_k24, min_k24]
    legend = ['max k24', 'min k24']

    if setting == 0:
        # plot simulations to check timecourse differences
        for i, j, k in zip(variables, color, legend):
            estimated_parameters = i[[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
                                      17, 18, 19, 20, 21, 22, 23, 24, 25, 30, 32]]
            estimated_species = i[[26, 27, 28, 29, 31]]
            parameters = shuffle.assign_estimated_parameters(estimated_parameters, prior_parameters)
            species = shuffle.assign_estimated_species(estimated_species, prior_species)
            parameters_all, species_all = calculate.dependent_values(species, parameters)
            sol = solve_ivp(simulatejakstat, [0, np.max(tspan_fine)], species_all.flatten(), method='BDF',
                            args=(parameters_all.flatten(),), t_eval=tspan_fine, rtol=1e-9, atol=1e-12)
            simulation = np.transpose(sol.y)
            SOCSindex1 = np.where(species_dictionary == 'SOCS1')[0][0]
            SOCSindex2 = np.where(species_dictionary == 'SOCS1PRLRJ2a')[0][0]
            SOCSindex3 = np.where(species_dictionary == 'SOCS1PRLRJ2aSHP2')[0][0]
            total_SOCS1 = simulation[:, SOCSindex1] + simulation[:, SOCSindex2] + simulation[:, SOCSindex3]
            normalized_SOCS1 = np.divide(total_SOCS1, np.max(total_SOCS1))
            plt.plot(np.divide(tspan_fine, 60), normalized_SOCS1, color=j)
            plt.ylabel('normalized total SOCS')
            plt.xlabel('minutes')
            plt.legend(['max k24', 'min k24'])
            plt.show()
            plt.savefig('posteriorapproximation/verifying_SOCS_timecourse_initializations.pdf')

    # generate initializations in correct format
    for i, k in zip(variables, legend):
        # simulate
        estimated_parameters = i[[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
                                  17, 18, 19, 20, 21, 22, 23, 24, 25, 30, 32]]
        estimated_species = i[[26, 27, 28, 29, 31]]
        parameters = shuffle.assign_estimated_parameters(estimated_parameters, prior_parameters)
        species = shuffle.assign_estimated_species(estimated_species, prior_species)
        parameters_all, species_all = calculate.dependent_values(species, parameters)
        sol = solve_ivp(simulatejakstat, [0, np.max(tspan_fine)], species_all.flatten(), method='BDF',
                        args=(parameters_all.flatten(),), t_eval=tspan_fine, rtol=1e-9, atol=1e-12)
        simulation_fine = np.transpose(sol.y)
        simulation = simulation_fine[indices, :]

        # calculate experimental predictions
        pStatA_norm_prop, pStatB_norm_prop, transloc_predA_prop, transloc_predB_prop, BcL_foldchange_prop, pJAK2_norm_prop \
            = calculate.experimental_predictions(simulation)
        predictions = np.concatenate((BcL_foldchange_prop, pJAK2_norm_prop, pStatA_norm_prop,
                                      pStatB_norm_prop, transloc_predA_prop, transloc_predB_prop), axis=None)
        # log likelihood
        log_likelihood = likelihood(predictions)

        # format and save initialization matrix
        initialization = np.array([parameters, species, log_likelihood, predictions, simulation_fine], dtype=object)
        np.save('posteriorapproximation/initialization_matrix_{}.npy'.format(k), initialization)


def metropolishastingsalgo(initialization, chainlength, prior_std, proposal_std, likelihood, setting):
    #initialize arrays
    parameters = initialization[0]
    species = initialization[1]
    likelihood_sc = initialization[2]
    exp_pred = initialization[3]
    # sim = initialization[4]
    if setting == 0:
        # initialize parameter, species, likelihood, simulation, acceptance rate lists to save
        biochemical_parameters = [parameters, parameters]
        biochemical_species = [species, species]
        likelihood_score = [likelihood_sc, likelihood_sc]
        exp_predictions = [exp_pred, exp_pred]
        # simulations = [sim, sim]
        acceptance_rate = []
        acceptance_count = []
    if setting == 1:
        # initialize parameter, species, likelihood, simulation, acceptance rate lists to save
        biochemical_parameters = list(parameters)
        biochemical_species = list(species)
        likelihood_score = list(likelihood_sc)
        exp_predictions = list(exp_pred)
        # simulations = list(sim)
        acceptance_rate = list(initialization[4])
        acceptance_count = []

    # initialize initial guess for parameter&species matrix
    fixed_params = np.load('posteriorapproximation/prior_parameters.npy')
    fixed_species = np.load('posteriorapproximation/prior_species.npy')

    #initialize tspan for ODE simulation
    tspan_course = np.load('posteriorapproximation/tspan_course_all_sec.npy')
    tspan_fine = np.load('posteriorapproximation/tspan_fine_all_sec.npy')
    indices = []
    for i in tspan_course:
        index = np.argwhere(tspan_fine == i)
        indices.append(index)
    indices = np.array(indices).flatten()

    for i in np.arange(0, chainlength):
        print(i)
        # generate proposed parameterization
        estimated_parameters = shuffle.index_estimated_parameters(biochemical_parameters[-1])
        estimated_species = shuffle.index_estimated_species(biochemical_species[-1])
        sampled_estimated_parameters = generate_lognormal(estimated_parameters, proposal_std)
        sampled_estimated_species = generate_lognormal(estimated_species, proposal_std)
        proposed_parameters = shuffle.assign_estimated_parameters(sampled_estimated_parameters, fixed_params)
        proposed_species = shuffle.assign_estimated_species(sampled_estimated_species, fixed_species)
        proposed_parameters, proposed_species = calculate.dependent_values(proposed_species, proposed_parameters)

        # simulate with proposed parameters
        sol = solve_ivp(simulatejakstat, [0, np.max(tspan_course)], proposed_species.flatten(), method='BDF',
                        args=(proposed_parameters.flatten(),), t_eval=tspan_course, rtol=1e-9, atol=1e-12)
        proposed_simulation = np.transpose(sol.y)

        # calculate experimental predictions'
        pStatA_norm_prop, pStatB_norm_prop, transloc_predA_prop, transloc_predB_prop, BcL_foldchange_prop, pJAK2_norm_prop\
            = calculate.experimental_predictions(proposed_simulation)
        proposed_predictions = np.concatenate((BcL_foldchange_prop, pJAK2_norm_prop, pStatA_norm_prop,
                                               pStatB_norm_prop, transloc_predA_prop, transloc_predB_prop), axis=None)

        # previous experimental predictions
        previous_predictions = exp_predictions[-1]

        # calculate acceptance ratio
        # simplified log prior probabilities. lognormal shape chosen based on distribution of kcat in Bar Even 2011
        sampled_estimated_all_variables = np.append(np.array(sampled_estimated_parameters), np.array(sampled_estimated_species))
        estimated_all_variables = np.append(np.array(estimated_parameters), np.array(estimated_species))
        initial_parameters = shuffle.index_estimated_parameters(fixed_params)
        initial_species = shuffle.index_estimated_species(fixed_species)
        initial_estimated_all_variables = np.append(np.array(initial_parameters), np.array(initial_species))
        log_prior_probability_proposed = evaluate_lognormal_function(sampled_estimated_all_variables, initial_estimated_all_variables, prior_std)
        log_prior_probability_previous = evaluate_lognormal_function(estimated_all_variables, initial_estimated_all_variables, prior_std)
        # simplified log proposal probabilities; lognormal shape chosen based on support over posterior'
        log_probability_proposed_params_given_previous = evaluate_lognormal_function(sampled_estimated_all_variables, estimated_all_variables, proposal_std)
        log_probability_previous_params_given_proposed = evaluate_lognormal_function(estimated_all_variables, sampled_estimated_all_variables, proposal_std)
        # log likelihood
        log_likelihood_proposed = likelihood(proposed_predictions)
        print('log likelihood', log_likelihood_proposed)
        log_likelihood_previous = likelihood(previous_predictions)
        # log acceptance probability
        log_acceptance_probability = log_prior_probability_proposed+log_likelihood_proposed+\
                                     log_probability_previous_params_given_proposed - \
                                     log_prior_probability_previous - log_likelihood_previous - \
                                     log_probability_proposed_params_given_previous

        # compare acceptance ratio to random sample
        u = np.random.random()
        u_log = math.log(u)

        # update posterior sample, save metrics and simulations
        if u_log <= min(math.log(1), log_acceptance_probability):
            biochemical_parameters.append(proposed_parameters)
            biochemical_species.append(proposed_species)
            likelihood_score.append(log_likelihood_proposed)
            exp_predictions.append(proposed_predictions)
            # simulations.append(proposed_simulation_fine)
            acceptance_count.append(1)
        else:
            biochemical_parameters.append(biochemical_parameters[-1])
            biochemical_species.append(biochemical_species[-1])
            likelihood_score.append(log_likelihood_previous)
            exp_predictions.append(previous_predictions)
            # simulations.append(simulations[-1])
            acceptance_count.append(0)

        # Acceptance Rate Display
        if len(acceptance_count) >= 100 and i % 50 == 0:
            acceptance_rate_i = sum(acceptance_count[i-100:i])/100
            print('acceptance rate, last 100 runs:', acceptance_rate_i)
            acceptance_rate.append(acceptance_rate_i)

    return likelihood_score, biochemical_species, biochemical_parameters, acceptance_rate, exp_predictions, # simulations
