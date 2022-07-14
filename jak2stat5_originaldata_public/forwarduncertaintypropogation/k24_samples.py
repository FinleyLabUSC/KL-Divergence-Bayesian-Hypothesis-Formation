import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def visualize_k24_samples(k24_values, max_k24_samples, min_k24_samples, path):
    plt.figure()
    sns.histplot(np.log(k24_values), color='black')
    sns.histplot(np.log(max_k24_samples), color='orange')
    sns.histplot(np.log(min_k24_samples), color='darkturquoise')
    plt.xlabel('ln(k24)', fontsize=14)
    plt.xlim([-15,0])
    plt.ylabel('')
    # plt.legend(['all k24 samples', 'larger k24 samples', 'smaller k24 samples'])
    plt.savefig('{}/k24_samples.png'.format(path), dpi=300)

