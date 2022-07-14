from shape_classification import shape_classification
from example_socs_timecourse import example_socs_timecourse_4curves
from shape_classification import adjustable_shape_classification
import numpy as np

socssimulation = np.load('generateinputs/fine_grain_SOCS_simulation.npy')

#shape_classification function does not allow adjustable inputs, but saves shape classification for supplementary figure s4
# adds text to squares for figure 5A, and saves line plot of all classified simulations
shape_classification(socssimulation)

#generate supplementary figure s4
example_socs_timecourse_4curves(np.load('shapeclassification/h10_simulations.npy'),
                                np.load('shapeclassification/h01_simulations.npy'),
                                np.load('shapeclassification/h00_simulations.npy'),
                                np.load('shapeclassification/h11_simulations.npy'),
                                'shapeclassification')

#generate supplementary figure S5A
#adjustable_shape_classification takes three inputs: simulations to classify, short-term window, and sustained/transient threshold
adjustable_shape_classification(socssimulation, 5, 0.5)

