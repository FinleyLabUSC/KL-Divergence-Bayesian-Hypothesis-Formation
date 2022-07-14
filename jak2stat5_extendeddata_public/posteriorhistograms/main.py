from posteriorhistogramplot import histogramsposterior_logscale
from posteriorhistogramplot import histogramspriorposterior_logscale
import numpy as np

posterior_samples = np.load('generateinputs/all_chains_baseline and timepoint_proposal 0.1_estimatedvariables.npy')
path = 'posteriorhistograms'

#visualize posterior samples
histogramsposterior_logscale(posterior_samples, path)

#visualize posterior samples compared to prior
histogramspriorposterior_logscale(posterior_samples, np.sqrt(2), path)


