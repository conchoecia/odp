"""
These are the plotting functions used by many ODP programs
"""

import matplotlib

def format_matplotlib():
    """format the fonts and print options for the plots"""
    font = {'family' : 'sans-serif',
            'sans-serif' : 'DejaVu Sans', # removing this after finding that many users don't have Helvetica installed. :( https://github.com/conchoecia/odp/issues/34
            'weight' : 'normal',
            'size'   : 12}

    matplotlib.rc('font', **font)

    grid = {"color": ".95", "linestyle": "-"}
    # grid style
    matplotlib.rc('grid', **grid)

    # Preserve the vertical order of embedded images:
    matplotlib.rcParams['image.composite_image'] = False
    # text as font in pdf
    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42

def plot_decay(datastruct, outpath, outtsv):
    """
    This plots the decay of an ALG between number of genes in the main chromosome,
    and the number of genes in smaller chromosomes

    Parameters:
      - 0th datastruct - the data structure described below with the data to plot
      - 1st outpath    - the path to save the plot to
      - 2nd outtsv     - the path to save the processed data to

    The input is this datastructure:
      ALG dataframe
       - key: the ALG name
         - 0th element is a list of genes that are on orthologous chroms
         - 1st element is a list of genes that are not on orthologous chroms
           - The key for these dicts is the scaffold name
           - The value for both of these dicts is the count of genes on that scaffold in that category
    """
    import matplotlib.pyplot as plt
    # convert the datastruct to a simple dataframe
    ALGs = list(datastruct.keys())
    conserved = [sum(datastruct[x][0].values()) for x in ALGs]
    scattered = [sum(datastruct[x][1].values()) for x in ALGs]
    total     = [conserved[i] + scattered[i] for i in range(len(ALGs))]
    # turn those columns into a dataframe
    df = pd.DataFrame({"ALG":ALGs, "conserved":conserved, "scattered":scattered, "total":total})
    df = df.sort_values(by="total", ascending=False)
    df = df.reset_index(drop=True)
    df.to_csv(outtsv, sep="\t", index=False)

    fig, ax1 = plt.subplots()

    color = 'tab:red'
    ax1.set_xlabel('ALG size')
    ax1.set_ylabel('Distribution size', color=color)
    ax1.tick_params(axis='y', labelcolor=color)

    for index, row in df.iterrows():
        x1 = row["total"]
        x2 = row["total"]
        y1 = 0
        y2 = row["total"]
        ax1.plot([x1,x2],[y1,y2],'k-')
    ax1.plot(df["total"], df["total"],'ro')
    ax1.plot(df["total"], df["scattered"],'bo')

    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    ax2.set_ylim([0, 100])

    color = 'tab:blue'
    ax2.set_ylabel('percent conserved on ALGs', color=color)  # we already handled the x-label with ax1
    ax2.plot(df["total"], 100*(df["conserved"]/df["total"]), "b-")
    ax2.tick_params(axis='y', labelcolor=color)

    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.savefig(outpath)
