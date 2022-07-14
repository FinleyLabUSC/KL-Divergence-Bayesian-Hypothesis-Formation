import numpy as np
import matplotlib.pyplot as plt

def SOCS_k24_forward_uncertainty_propogation(max_k24, min_k24, path, percentile_lower, percentile_upper, fontsize):
    # initialize tspan for plot
    tspan_course = np.load('{}/tspan_course_all_sec.npy'.format(path))
    tspan_fine = np.load('{}/tspan_fine_all_sec.npy'.format(path))
    indices = []
    for i in tspan_course:
        index = np.argwhere(tspan_fine == i)
        indices.append(index)
    indices = np.array(indices).flatten()

    plt.figure(figsize=(11, 7))
    plt.title('SOCS relative change {}th to {}th percentile of simulations'.format(percentile_lower*100, percentile_upper*100))
    median_upper = np.median(np.array(max_k24), axis=0)
    median_lower = np.median(np.array(min_k24), axis=0)
    lower_bound_larger_k24 = np.array(np.quantile(max_k24, percentile_lower, axis=0), dtype=float)
    upper_bound_larger_k24 = np.array(np.quantile(max_k24, percentile_upper, axis=0), dtype=float)
    lower_bound_smaller_k24 = np.array(np.quantile(min_k24, percentile_lower, axis=0), dtype=float)
    upper_bound_smaller_k24 = np.array(np.quantile(min_k24, percentile_upper, axis=0), dtype=float)
    plt.subplot(2, 1, 1)
    plt.fill_between(tspan_fine, lower_bound_larger_k24, upper_bound_larger_k24, alpha=0.1, color='orange', label='faster degradation regime')
    plt.plot(tspan_fine, median_upper, color='orange', linestyle='dotted', linewidth=3)
    plt.xticks(ticks=tspan_course[[0, 8, len(tspan_course)-1]], labels=['','',''], fontsize=fontsize-1)
    #plt.xticks(ticks=tspan_course[4:len(tspan_course)], labels=['','','','','','','','','','','',''], rotation=90)
    # plt.xlabel('time, hours')
    plt.ylabel('Relative [SOCS]', fontsize=fontsize)
    #plt.legend(loc='lower center', fontsize=fontsize)
    plt.subplot(2, 1, 2)
    plt.fill_between(tspan_fine, lower_bound_smaller_k24, upper_bound_smaller_k24, alpha=0.1, color='darkturquoise', label='slower degradation regime')
    plt.plot(tspan_fine, median_lower, color='darkturquoise', linestyle='dotted', linewidth=3)
    #plt.legend(loc='lower center', fontsize=fontsize)
    plt.xticks(ticks=tspan_course[[0, 8, len(tspan_course) - 1]], labels=np.around(tspan_course[[0, 8, len(tspan_course) - 1]] / 3600, decimals=0), fontsize=fontsize-1)
    plt.xlabel('Time, hours', fontsize=fontsize)
    plt.ylabel('Relative [SOCS]', fontsize=fontsize)
    plt.savefig('{}/k24_sorted_SOCS_timecourse_{}th to {}th percentile.png'.format(path, percentile_lower*100, percentile_upper*100), dpi=300)
