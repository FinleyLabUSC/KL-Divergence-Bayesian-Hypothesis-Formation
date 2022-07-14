import numpy as np
'#Calculate dependent parameters and species after sampling proposal distribution'

def dependent_values(proposed_species, proposed_parameters):

    k1 = proposed_parameters[9]*proposed_species[1]
    # k1 = kdeg*RJ
    # synthesis rate of RJ
    k35 = proposed_parameters[59]*proposed_species[55]
    # k35 = k34*BcL
    s5ac = proposed_parameters[53]/(1+proposed_parameters[51]*proposed_parameters[50])
    # s5ac = totalSTAT/(1 + ncratioA*Vratio)
    # STAT5A in cytosol
    s5bc = proposed_parameters[53]/(1+proposed_parameters[52]*proposed_parameters[50])
    # s5bc = #totalSTAT/(1 + ncratioB*Vratio)
    # STAT5B in cytosol
    s5an = (proposed_parameters[53]-s5ac)/proposed_parameters[50]
    # s5an = (totalSTAT - S5Ac)/Vratio
    # STAT5A in nucleus
    s5bn = (proposed_parameters[53]-s5bc)/proposed_parameters[50]
    # s5bn = #(totalSTAT - S5Bc)/Vratio
    # STAT5B in nucleus
    k17ina = (s5an/s5ac)*proposed_parameters[50]*proposed_parameters[32]
    # k17ina = S5An/S5Ac*(Vratio)*k17outA
    k17inb = (s5bn/s5bc)*proposed_parameters[50]*proposed_parameters[34]
    # k17inb = S5Bn/S5Bc*(Vratio)*k17outB

    proposed_parameters[0] = k1
    proposed_parameters[60] = k35
    proposed_parameters[31] = k17ina
    proposed_parameters[33] = k17inb
    proposed_species[2] = s5ac
    proposed_species[3] = s5bc
    proposed_species[34] = s5an
    proposed_species[35] = s5bn

    return proposed_parameters,proposed_species

def experimental_predictions(simulation):
    '# Function to Calculate species of interest from ODE simulation. ODE must be run with particular tspan'
    # tspan should have length 16, and have 1 step = 1 sec

    predConc = simulation

    # pStat timecourse
    total_pStatA = predConc[:, 12] +2. * predConc[:, 14] + predConc[:, 15] + predConc[:, 18] + 2*predConc[:, 20] + \
               predConc[:, 22] + predConc[:, 23] + predConc[:, 26] + 2 * predConc[:, 27] + predConc[:, 28] + \
               predConc[:, 30] + predConc[:, 32] + 2 * predConc[:, 36] + predConc[:, 37] + predConc[:, 39] + \
               predConc[:, 42]
    total_pStatB = predConc[:, 13] +2. * predConc[:, 16] + predConc[:, 15] + predConc[:, 19] + 2 * predConc[:, 21] + \
               predConc[:, 22] + predConc[:, 24] + predConc[:, 25] + 2 * predConc[:, 29] + predConc[:, 28] + \
               predConc[:, 31] + predConc[:, 33] + 2 * predConc[:, 38] + predConc[:, 37] + predConc[:, 40] + \
               predConc[:, 41]

    pStatA_norm = [total_pStatA[0], total_pStatA[2], total_pStatA[4], total_pStatA[5], total_pStatA[7], total_pStatA[9],
               total_pStatA[11], total_pStatA[13]] / total_pStatA[4]
    pStatB_norm = [total_pStatB[0], total_pStatB[2], total_pStatB[4], total_pStatB[5], total_pStatB[7], total_pStatB[9],
               total_pStatB[11], total_pStatB[13]] / total_pStatB[4]

    # Stat translocation
    Stat_cytoA = predConc[:, 2] + predConc[:, 10] + predConc[:, 12] + 2 * predConc[:, 14] + predConc[:, 15] + \
                 predConc[:, 18] + 2 * predConc[:, 20] + predConc[:, 22] + 2 * predConc[:, 23] + predConc[:, 25] + \
                 predConc[:, 26] + predConc[:, 47]
    Stat_cytoB = predConc[:, 3] + predConc[:, 11] + predConc[:, 13] + 2 * predConc[:, 16] + predConc[:, 15] + \
                 predConc[:, 19] + 2 * predConc[:, 21] + predConc[:, 22] + 2 * predConc[:, 24] + predConc[:, 25] + \
                 predConc[:, 26] + predConc[:, 48]

    Stat_nucleusA = 2 * predConc[:, 27] + predConc[:, 28] + predConc[:, 30] + predConc[:, 32] + predConc[:, 34] + \
                    2 * predConc[:, 36] + predConc[:, 37] + 2 * predConc[:, 39] + predConc[:, 41] + predConc[:, 42]
    Stat_nucleusB = 2 * predConc[:, 29] + predConc[:, 28] + predConc[:, 31] + predConc[:, 33] + predConc[:, 35] + \
                    2 * predConc[:, 38] + predConc[:, 37] + 2 * predConc[:, 40] + predConc[:, 41] + predConc[:, 42]

    nucleus_cyto_ratioA = Stat_nucleusA / Stat_cytoA
    nucleus_cyto_ratioB = Stat_nucleusB / Stat_cytoB

    transloc_predA = [nucleus_cyto_ratioA[4], nucleus_cyto_ratioA[5], nucleus_cyto_ratioA[6], nucleus_cyto_ratioA[8],
                      nucleus_cyto_ratioA[11]]
    transloc_predB = [nucleus_cyto_ratioB[1], nucleus_cyto_ratioB[3], nucleus_cyto_ratioB[4], nucleus_cyto_ratioB[5],
                  nucleus_cyto_ratioB[7], nucleus_cyto_ratioB[8], nucleus_cyto_ratioB[9], nucleus_cyto_ratioB[10],
                  nucleus_cyto_ratioB[11]]

    # Bcl - xL foldchange
    BcL_foldchange = [predConc[0, 55], predConc[7, 55], predConc[9, 55], predConc[11, 55], predConc[12, 55], predConc[13, 55],
                         predConc[14, 55], predConc[15, 55]] / predConc[0, 55]

    # pJAK2
    total_pJAK2 = predConc[:, 9] + predConc[:, 10] + predConc[:, 11] + predConc[:, 17] + predConc[:, 46] + \
                  predConc[:, 47] + predConc[:, 48] + predConc[:, 49];
    pJAK2_norm = [total_pJAK2[0], total_pJAK2[2], total_pJAK2[4], total_pJAK2[5], total_pJAK2[7], total_pJAK2[9],
                  total_pJAK2[11], total_pJAK2[13]] / total_pJAK2[2]

    return pStatA_norm, pStatB_norm, transloc_predA, transloc_predB, BcL_foldchange, pJAK2_norm

def experimental_predictions_fine_tspan(predConc):
    #tspan is sampled every 20 seconds

    # pStat timecourse
    total_pStatA = predConc[:, 12] +2. * predConc[:, 14] + predConc[:, 15] + predConc[:, 18] + 2*predConc[:, 20] + \
               predConc[:, 22] + predConc[:, 23] + predConc[:, 26] + 2 * predConc[:, 27] + predConc[:, 28] + \
               predConc[:, 30] + predConc[:, 32] + 2 * predConc[:, 36] + predConc[:, 37] + predConc[:, 39] + \
               predConc[:, 42]
    total_pStatB = predConc[:, 13] +2. * predConc[:, 16] + predConc[:, 15] + predConc[:, 19] + 2 * predConc[:, 21] + \
               predConc[:, 22] + predConc[:, 24] + predConc[:, 25] + 2 * predConc[:, 29] + predConc[:, 28] + \
               predConc[:, 31] + predConc[:, 33] + 2 * predConc[:, 38] + predConc[:, 37] + predConc[:, 40] + \
               predConc[:, 41]

    pStatA_norm = np.divide(total_pStatA, total_pStatA[90])
    pStatB_norm = np.divide(total_pStatB,total_pStatB[90])

    # Stat translocation
    Stat_cytoA = predConc[:, 2] + predConc[:, 10] + predConc[:, 12] + 2 * predConc[:, 14] + predConc[:, 15] + \
                 predConc[:, 18] + 2 * predConc[:, 20] + predConc[:, 22] + 2 * predConc[:, 23] + predConc[:, 25] + \
                 predConc[:, 26] + predConc[:, 47]
    Stat_cytoB = predConc[:, 3] + predConc[:, 11] + predConc[:, 13] + 2 * predConc[:, 16] + predConc[:, 15] + \
                 predConc[:, 19] + 2 * predConc[:, 21] + predConc[:, 22] + 2 * predConc[:, 24] + predConc[:, 25] + \
                 predConc[:, 26] + predConc[:, 48]

    Stat_nucleusA = 2 * predConc[:, 27] + predConc[:, 28] + predConc[:, 30] + predConc[:, 32] + predConc[:, 34] + \
                    2 * predConc[:, 36] + predConc[:, 37] + 2 * predConc[:, 39] + predConc[:, 41] + predConc[:, 42]
    Stat_nucleusB = 2 * predConc[:, 29] + predConc[:, 28] + predConc[:, 31] + predConc[:, 33] + predConc[:, 35] + \
                    2 * predConc[:, 38] + predConc[:, 37] + 2 * predConc[:, 40] + predConc[:, 41] + predConc[:, 42]

    transloc_predA = Stat_nucleusA / Stat_cytoA
    transloc_predB = Stat_nucleusB / Stat_cytoB

    # Bcl - xL foldchange
    BcL_foldchange = np.divide(predConc[:,55], predConc[0,55])

    # pJAK2
    total_pJAK2 = predConc[:, 9] + predConc[:, 10] + predConc[:, 11] + predConc[:, 17] + predConc[:, 46] + \
                  predConc[:, 47] + predConc[:, 48] + predConc[:, 49];
    pJAK2_norm = np.divide(total_pJAK2, total_pJAK2[30])

    return pStatA_norm, pStatB_norm, transloc_predA, transloc_predB, BcL_foldchange, pJAK2_norm
