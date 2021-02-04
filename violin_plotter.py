"""
=========================
Violin plot customization
=========================

#### THE CODE BELOW CRIBS HEAVILY FROM THIS SOURCE: 
https://matplotlib.org/3.3.1/gallery/statistics/customized_violin.html#sphx-glr-gallery-statistics-customized-violin-py

This example demonstrates how to fully customize violin plots.
The first plot shows the default style by providing only
the data. The second plot first limits what matplotlib draws
with additional kwargs. Then a simplified representation of
a box plot is drawn on top. Lastly, the styles of the artists
of the violins are modified.

For more information on violin plots, the scikit-learn docs have a great
section: https://scikit-learn.org/stable/modules/density.html
"""

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm 

plt.rcParams["font.family"] = 'serif'


def create_violin(data, data_labels=None, x_label=None, y_label=None, plot_title=None, colormap='viridis', autoshow=False):
    #### data should be a list of lists!
    #### autoshow=True means the plot will show up as soon as you call the function. Leave as FALSE to make modifications in a script.

    def adjacent_values(vals, q1, q3):
        upper_adjacent_value = q3 + (q3 - q1) * 1.5
        upper_adjacent_value = np.clip(upper_adjacent_value, q3, vals[-1])

        lower_adjacent_value = q1 - (q3 - q1) * 1.5
        lower_adjacent_value = np.clip(lower_adjacent_value, vals[0], q1)
        return lower_adjacent_value, upper_adjacent_value


    def set_axis_style(ax, labels=None, x_label=x_label, y_label=y_label):
        ax.get_xaxis().set_tick_params(direction='out')
        ax.xaxis.set_ticks_position('bottom')
        if type(labels) != type(None):
            ax.set_xticks(np.arange(1, len(labels) + 1))
            ax.set_xticklabels(labels)
            ax.set_xlim(0.25, len(labels) + 0.75)
        try:
            ax.set_xlabel(x_label)
        except:
            pass

        try:
            ax.set_ylabel(y_label)
        except:
            pass


    # create test data
    #np.random.seed(19680801)
    #data = [sorted(np.random.normal(0, std, 100)) for std in range(1, 5)]

    fig, ax2 = plt.subplots(figsize=(9, 4))

    #ax1.set_title(plot_title)
    #ax1.set_ylabel(y_label)
    #ax1.violinplot(data)

    ax2.set_title(plot_title)
    parts = ax2.violinplot(
            data, showmeans=False, showmedians=False,
            showextrema=False)

    #colors = cm.inferno(np.linspace(0,1,np.array(data).shape[0])) 
    ncolors = np.array(data).shape[0]
    colors = cm.get_cmap(colormap)(np.linspace(0,1,ncolors))


    for npc,pc in enumerate(parts['bodies']):
        #pc.set_facecolor('#D43F3A')
        pc.set_facecolor(colors[npc])
        pc.set_edgecolor('black')
        pc.set_alpha(1)

    #quartile1, medians, quartile3 = np.nanpercentile(data, [25, 50, 75], axis=1)
    quartile1, medians, quartile3 = [], [], []
    for i in np.arange(0,ncolors,1):
        q1,med,q3 = np.nanpercentile(np.array(data[i]), [25,50,70])
        quartile1.append(q1)
        medians.append(med)
        quartile3.append(q3)
    quartile1, medians, quartile3 = np.array(quartile1), np.array(medians), np.array(quartile3)

    whiskers = np.array([
        adjacent_values(sorted_array, q1, q3)
        for sorted_array, q1, q3 in zip(data, quartile1, quartile3)])
    whiskers_min, whiskers_max = whiskers[:, 0], whiskers[:, 1]

    inds = np.arange(1, len(medians) + 1)
    ax2.scatter(inds, medians, marker='o', color='white', s=30, zorder=3)
    ax2.vlines(inds, quartile1, quartile3, color='k', linestyle='-', lw=5)
    ax2.vlines(inds, whiskers_min, whiskers_max, color='k', linestyle='-', lw=1)

    # set style for the axes
    #labels = ['A', 'B', 'C', 'D']
    labels = data_labels
    #for ax in [ax1, ax2]:
    #    set_axis_style(ax, labels)
    set_axis_style(ax2, labels)

    plt.subplots_adjust(left=0.09, bottom=0.15, right=0.92, top=0.88, wspace=0.05, hspace=0.05)
    if autoshow == True:
        plt.show()
