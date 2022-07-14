import numpy as np

def index_estimated_parameters(parameter_array):

    parameter_array = np.float64(parameter_array)
    '#ensure array is int64 type'

    mult8B = parameter_array[13]/parameter_array[11]
    '#mult8B = k8B/k8A'
    mult8AB = parameter_array[15]/parameter_array[11]
    '#mult8AB = k8AB/k8A'
    mult14B = parameter_array[26]/parameter_array[25]
    '#mult14B = k14B/k14A'
    mult14AB = parameter_array[27]/parameter_array[25]
    '#mult14B = k14AB/k14A'
    '#mult14AB = parameter_array[27]/parameter_array[11]'
    '#mult14AB = k14Ab/k8A'
    mult17B = parameter_array[34]/parameter_array[32]
    '#mult17B = k17Bout/k17Aout'

    estimated_parameters = parameter_array[[3,6,10,11,17,20,24,25,28,30,32,37,39,41,42,43,44,47,48,53,54,59,5]]
    estimated_parameters = np.insert(estimated_parameters,4,[mult8B])
    estimated_parameters = np.insert(estimated_parameters, 5, [mult8AB])
    estimated_parameters = np.insert(estimated_parameters, 10, [mult14B])
    estimated_parameters = np.insert(estimated_parameters,11,[mult14AB])
    estimated_parameters = np.insert(estimated_parameters, 15, [mult17B])
    return estimated_parameters

def assign_estimated_parameters(sampled_estimated_parameters, parameters):

    parameter_array = np.arange(0,61)
    parameter_array = np.float64(parameter_array)
    parameters = np.float64(parameters)
    '#ensure array is int64 type'
    sampled_estimated_parameters = np.float64(sampled_estimated_parameters)
    '#ensure array is int64 type'
    parameter_array[13] = sampled_estimated_parameters[4] * sampled_estimated_parameters[3]
    'k8B = mult8B*k8A'
    parameter_array[15] = sampled_estimated_parameters[5] * sampled_estimated_parameters[3]
    'k8AB = mult8AB*k8A'
    parameter_array[26] = sampled_estimated_parameters[10] * sampled_estimated_parameters[9]
    'k14B = mult14B*k14A'
    parameter_array[27] = sampled_estimated_parameters[11] * sampled_estimated_parameters[9]
    'k14AB = mult14AB*k14A'
    parameter_array[34] = sampled_estimated_parameters[15] * sampled_estimated_parameters[14]
    'k17outB = mult17B*k17outA'

    for i,j in zip([3,6,10,11,17,20,24,25,28,30,32,37,39,41,42,43,44,47,48,53,54,59,5],
                   [0,1,2,3,6,7,8,9,12,13,14,16,17,18,19,20,21,22,23,24,25,26,27]):
        parameter_array[i] = sampled_estimated_parameters[j]
    for i in [0,1,2,4,7,8,9,12,14,16,18,19,21,22,23,29,31,33,35,36,38,40,45,46,49,
                    50,51,52,55,56,57,58,60]:
        parameter_array[i] = parameters[i]
    return parameter_array

def index_estimated_species(species_array):
    species_array = np.float64(species_array)
    '#ensure array is int64 type'
    return species_array[[1,4,5,6,55]]

def assign_estimated_species(sampled_estimated_species,species):
    species_array = np.arange(0,56)
    species_array = np.float64(species_array)
    species = np.float64(species)
    '#ensure array is int64 type'
    sampled_estimated_species = np.float64(sampled_estimated_species)
    '#ensure array is int64 type'

    for i,j in zip([1,4,5,6,55],np.arange(0,5,1)):
        species_array[i]=sampled_estimated_species[j]
    for i in [0,2,3,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,
              41,42,43,44,45,46,47,48,49,50,51,52,53,54]:
        species_array[i]=species[i]
    return species_array

def index_estimated_parameters_and_species(parameter_matrix,species_matrix):
    estimated_parameters = parameter_matrix[:,[3,6,10,11,17,20,24,25,28,30,32,37,39,41,42,43,44,47,48,53,54,59,5]]
    mult8B = np.divide(parameter_matrix[:,13],parameter_matrix[:,11])
    '#mult8B = k8B/k8A'
    mult8AB = np.divide(parameter_matrix[:,15],parameter_matrix[:,11])
    '#mult8AB = k8AB/k8A'
    mult14B = np.divide(parameter_matrix[:,26],parameter_matrix[:,25])
    '#mult14B = k14B/k14A'
    mult14AB = np.divide(parameter_matrix[:,27],parameter_matrix[:,25])
    '#mult14B = k14AB/k14A'
    mult17B = np.divide(parameter_matrix[:,34],parameter_matrix[:,32])
    '#mult17B = k17Bout/k17Aout'
    estimated_parameters = np.insert(estimated_parameters, 4,[mult8B],axis=1)
    estimated_parameters = np.insert(estimated_parameters, 5, [mult8AB],axis=1)
    estimated_parameters = np.insert(estimated_parameters, 10, [mult14B],axis=1)
    estimated_parameters = np.insert(estimated_parameters, 11,[mult14AB],axis=1)
    estimated_parameters = np.insert(estimated_parameters, 15, [mult17B],axis=1)
    estimated_species = species_matrix[:, [1, 4, 5, 6, 55]]
    estimated_variables = np.concatenate((estimated_parameters,estimated_species),axis=1)
    return estimated_variables

def isolate_estimated_variables(posteriorparameters, posteriorspecies, path_name):
    '#import chains'
    parameterchain1 = posteriorparameters
    specieschain1 = posteriorspecies

    # take off burn-in
    parameterchain1 = parameterchain1[1000:-1, :]
    specieschain1 = specieschain1[1000:-1, :]

    # index estimated variables
    parameterchain1est = []
    specieschain1est = []

    for i in np.arange(0, len(parameterchain1)):
        estimated_parameters1 = index_estimated_parameters(parameterchain1[i, :])
        estimated_species1 = index_estimated_species(specieschain1[i, :])
        parameterchain1est.append(estimated_parameters1)
        specieschain1est.append(estimated_species1)

    parameterchain1est = np.array(parameterchain1est)
    specieschain1est = np.array(specieschain1est)

    '#combine estimated variables into one matrix for back-compatibility'
    estimated_parameters1 = np.concatenate((parameterchain1est[:, 0:26], specieschain1est[:, 0:4]), axis=1)
    estimated_parameters1 = np.concatenate((estimated_parameters1, np.transpose([parameterchain1est[:, 26], ])), axis=1)
    estimated_parameters1 = np.concatenate((estimated_parameters1, np.transpose([specieschain1est[:, 4]])), axis=1)
    estimated_parameters1 = np.concatenate((estimated_parameters1, np.transpose([parameterchain1est[:, 27]])), axis=1)

    np.save('{}'.format(path_name), estimated_parameters1)

def isolate_opposite_sides_k24(estimated_variables, estimated_variables_dictionary):
    k24_index = np.where(estimated_variables_dictionary == 'k24')[0][0]
    estimated_variables_sorted = estimated_variables[estimated_variables[:, k24_index].argsort()]
    max_k24 = estimated_variables_sorted[-1, :]
    min_k24 = estimated_variables_sorted[0, :]
    return max_k24,  min_k24




