import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe

def example_socs_timecourse(max_k24, path):
    # initialize tspan for plot
    tspan_course = np.load('{}/tspan_course_all_sec.npy'.format(path))
    tspan_fine = np.load('{}/tspan_fine_all_sec.npy'.format(path))
    indices = []
    for i in tspan_course:
        index = np.argwhere(tspan_fine == i)
        indices.append(index)
    indices = np.array(indices).flatten()
    plt.figure(figsize=(11, 7))
    median_upper = np.median(np.array(max_k24), axis=0)
    plt.plot(tspan_fine, median_upper, color='black', linewidth=3)
    plt.xticks(ticks=tspan_course[[0, 10, len(tspan_course)-1]], labels=np.around(tspan_course[[0, 10, len(tspan_course)-1]] / 3600, decimals=0),
               fontsize=14)
    plt.yticks([0, 0.5, 1], fontsize=14)
    plt.xlabel('time, hours', fontsize=14)
    plt.ylabel('relative [SOCS]', fontsize=14)
    plt.savefig('{}/example_SOCS_timecourse.png'.format(path), dpi=300)

def example_socs_timecourse_4curves(h10,h01,h00,h11, path):
    # initialize tspan for plot
    tspan_course = np.load('{}/tspan_course_all_sec.npy'.format(path))
    tspan_fine = np.load('{}/tspan_fine_all_sec.npy'.format(path))
    indices = []
    for i in tspan_course:
        index = np.argwhere(tspan_fine == i)
        indices.append(index)
    indices = np.array(indices).flatten()
    plt.figure(figsize=(11, 7))
    color = ['darkorchid', 'silver', 'whitesmoke', 'dimgrey']
    shape_simulations = [h10, h01, h00, h11]
    plt.fill_betweenx([0,1], 0, tspan_course[10], color='dodgerblue', alpha=0.3)
    #plt.fill_betweenx([0, 1], tspan_fine[-5], tspan_course[len(tspan_course)-1], color='firebrick')
    plt.hlines(0.5, 0, tspan_fine[-1], linestyles='dashed', color='firebrick', linewidth=3)
    plt.vlines(tspan_fine[-1],0,1, linestyles='solid', color='firebrick', linewidth=3)
    for i in np.arange(0, 4):
        if i == 2:
            plt.plot(tspan_fine, np.median(shape_simulations[i], axis=1), color=color[i], linewidth=3,path_effects=[pe.Stroke(linewidth=5, foreground='k'), pe.Normal()])
        else:
            plt.plot(tspan_fine, np.median(shape_simulations[i], axis=1), color=color[i], linewidth=3)
    plt.xticks(ticks=tspan_course[[0, 10, len(tspan_course)-1]], labels=np.around(tspan_course[[0, 10, len(tspan_course)-1]] / 3600, decimals=0),
               fontsize=14)
    plt.yticks([0, 0.5, 1], fontsize=14)
    plt.xlabel('Time, hours', fontsize=14)
    plt.ylabel('Relative [SOCS]', fontsize=14)
    plt.savefig('{}/example_SOCS_timecourse_4curves.png'.format(path), dpi=300)
