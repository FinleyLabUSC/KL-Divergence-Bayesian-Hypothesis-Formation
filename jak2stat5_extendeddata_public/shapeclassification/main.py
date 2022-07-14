from shape_classification import shape_classification
from shape_classification import adjustable_shape_classification
import numpy as np


socssimulation = np.load('generateinputs/fine_grain_SOCS_simulation.npy')

#shape_classification does not allow adjustable inputs, but adds text to squares for figure 5B
# and saves line plot of all classified simulations
shape_classification(socssimulation)

#generate supplementary figure S5B
#adjustable_shape_classification takes three inputs: simulations to classify, short-term window, and sustained/transient threshold
adjustable_shape_classification(socssimulation, 5, 0.7)