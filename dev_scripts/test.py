#!/usr/bin/env python3
import pandas as pd
import sys

init = "random_sim_init_plottable.txt"
data = "random_sim_final_plottable_COW.txt.head"

import matplotlib.cm as cm
from matplotlib import pyplot as plt
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
initdf = pd.read_csv(init, index_col = None, sep = "\t")
df = pd.read_csv(data, index_col = None, sep = "\t")

import seaborn as sns
sns.set_theme(style="ticks")

possible_vars = list(set(["_".join(x.split("_")[1::]) for x in df.columns]))
plot_these = ["total_number_of_genes_in_FLGs",
              "total_number_of_groupings",
              "max_number_of_genes_in_grouping",
              "max_number_of_genes_in_single_FLG"]

plot_to_xlabel = {"total_number_of_genes_in_FLGs": "Genes in PI LGs",
                  "total_number_of_groupings": "Number of LG groupings",
                  "max_number_of_genes_in_grouping": "Genes in largest LG grouping",
                  "max_number_of_genes_in_single_FLG": "Genes in largest LG"
                 }


print(plot_these)
print(possible_vars)
for x in plot_these:
    if x not in possible_vars:
        raise IOError("{} was not in the input df".format(x))

for thisthing in plot_these:
    fig, ax = plt.subplots(2, 1, sharey = True, sharex = True)
    fig.suptitle("{} randomized".format(thisthing), fontsize=8)
    fig.set_size_inches(5, 2.5)
    columns = ["{}_{}".format(x, thisthing) for x in ["hyp1", "hyp2"]]
    index = 0

    hyp1_init = initdf.iloc[0][columns[0]]
    hyp2_init = initdf.iloc[0][columns[1]]
    xmax = int(max([hyp1_init, hyp2_init,
                max(df[columns[0]]),
                max(df[columns[1]])
              ]) * 1.2)
    ymax = max([len(df[x]) for x in columns])
    for thiscol in columns:
        sns.histplot(df[thiscol], discrete = True,
                     label=thiscol, color = "#000000",
                          ax=ax[index])
        ax[index].set_title(thiscol, fontsize = 6)
        ax[index].ticklabel_format(axis='y', style='sci', scilimits=(0,0))
        ax[index].tick_params(axis='both', which='major', labelsize=8)
        ax[index].xaxis.label.set_size(8)
        ax[index].yaxis.label.set_size(8)
        sns.despine(offset=2, trim=False)
        index += 1

    # set the plotting limits
    #ax[0].set_xlim([-1,xmax])
    #ax[1].set_xlim([-1,xmax])
    ylim = ax[0].get_ylim()[-1]
    ax[0].axvline(hyp1_init, 0, (ylim*0.75)/ylim, color = "#c44e52")
    ax[1].axvline(hyp2_init, 0, (ylim*0.75)/ylim, color = "#c44e52")

    # set the xlabel
    ax[1].set_xlabel(plot_to_xlabel[thisthing])
    plt.tight_layout()
    plt.show()

    #fig = g.fig
    #pdf_pages.savefig(fig)
    #plt.close(g.fig)
