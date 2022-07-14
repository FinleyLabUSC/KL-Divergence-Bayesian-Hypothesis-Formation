from posteriorhistogramplot import histogramsposterior_logscale
from posteriorhistogramplot import histogramspriorposterior_logscale
from posteriorhistogramplot import histogramspriorposterior_logscale_top_kl_div
import numpy as np

posterior_samples = np.load('generateinputs/all_chains_baseline_proposal 0.1_estimatedvariables.npy')
path = 'posteriorhistograms'

#visualize posterior samples
histogramsposterior_logscale(posterior_samples, path)

#visualize posterior samples compared to prior
histogramspriorposterior_logscale(posterior_samples, np.sqrt(2), path)

#generate figure 3B
histogramspriorposterior_logscale_top_kl_div(posterior_samples, np.sqrt(2), path, 'black')

