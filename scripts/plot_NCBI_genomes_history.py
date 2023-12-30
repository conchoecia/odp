#!/usr/bin/env python

"""
This script plots the history of the availability of genome assemblies.

Makes separate plots for the availability of assemblies:
  1. Given my filtering criteria, specified in the rule filter_raw_genome_df in the snakefile GenDB_scrape_genomes_NCBI.snakefile
  2. The raw number of chromosome scale and annotated assemblies

The data table looks like this:
        coltype chrscale annotated  num_species  num_genomes assembly_release_date_cutoff
        cell    False     False         4745         4745                   2023-12-29
        cell    False      True          851          851                   2023-12-29
        cell     True     False         2163         2837                   2023-12-29
        cell     True      True          648          789                   2023-12-29
    marginal    False       all         5248         5596                   2023-12-29
    marginal     True       all         2287         3626                   2023-12-29
    marginal      all     False         6908         7582                   2023-12-29
    marginal      all      True         1499         1640                   2023-12-29
       total      all       all         7535         9222                   2023-12-29
"""

import argparse
from datetime import datetime
import matplotlib.pyplot as plt
from matplotlib import dates as mdates
from matplotlib.backends.backend_pdf import PdfPages
import os
import pandas as pd
import sys

def parse_args():
    """
    The only args we need are a path to the data table and a path to the output pdf.
    """
    parser = argparse.ArgumentParser(description='Plot the history of the availability of genome assemblies.')
    parser.add_argument('-i', '--input', help='Path to the data table.')
    parser.add_argument('-o', '--output', help='Path to the output pdf.')
    args = parser.parse_args()
    # make sure that the input file exists
    if not os.path.isfile(args.input):
        raise ValueError('The input file {} does not exist.'.format(args.input))
    return args

def lineplot_of_my_filtered_assemblies(df) -> plt.Figure:
    """
    Returns a matplotlib figure object.
    This is a two-panel top-bottom lineplot of the "cell" coltypes.
    For each plot there is also the absolute total from the "coltype == total" row.
      - The top panel is the number of genomes,
      - The bottom panel is the number of species.

    """
    # Set up a two-panel top-bottom line plot. Both plots will be rectangles and will share the x-axis.
    #  The x-axis will be the assembly release date cutoff, plotted as a scatterplot.
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(10, 7))
    # groupby the assembly release date cutoff, sort by date, descending
    gb = df.groupby('assembly_release_date_cutoff')

    xaxis_plotdate = []
    yaxis_top_genomes = {"chrscale_annotated": [],
                         "chrscale_unannotated": [],
                         "nonchrscale_annotated": [],
                         "nonchrscale_unannotated": [],
                         "total": []
                         }
    yaxis_bot_species = {"chrscale_annotated": [],
                         "chrscale_unannotated": [],
                         "nonchrscale_annotated": [],
                         "nonchrscale_unannotated": [],
                         "total": []
                         }
    for name, group in gb:
        # convert the name to a datetime object
        xaxis_plotdate.append(datetime.strptime(name, '%Y-%m-%d'))
        # add the value from this date to each of the yaxis lists
        yaxis_top_genomes["chrscale_annotated"].append(group.loc[(group['coltype'] == 'cell') & (group['chrscale'] == True) & (group['annotated'] == True), 'num_genomes'].values[0])
        yaxis_top_genomes["chrscale_unannotated"].append(group.loc[(group['coltype'] == 'cell') & (group['chrscale'] == True) & (group['annotated'] == False), 'num_genomes'].values[0])
        yaxis_top_genomes["nonchrscale_annotated"].append(group.loc[(group['coltype'] == 'cell') & (group['chrscale'] == False) & (group['annotated'] == True), 'num_genomes'].values[0])
        yaxis_top_genomes["nonchrscale_unannotated"].append(group.loc[(group['coltype'] == 'cell') & (group['chrscale'] == False) & (group['annotated'] == False), 'num_genomes'].values[0])
        yaxis_top_genomes["total"].append(group.loc[(group['coltype'] == 'total'), 'num_genomes'].values[0])

        yaxis_bot_species["chrscale_annotated"].append(group.loc[(group['coltype'] == 'cell') & (group['chrscale'] == True) & (group['annotated'] == True), 'num_species'].values[0])
        yaxis_bot_species["chrscale_unannotated"].append(group.loc[(group['coltype'] == 'cell') & (group['chrscale'] == True) & (group['annotated'] == False), 'num_species'].values[0])
        yaxis_bot_species["nonchrscale_annotated"].append(group.loc[(group['coltype'] == 'cell') & (group['chrscale'] == False) & (group['annotated'] == True), 'num_species'].values[0])
        yaxis_bot_species["nonchrscale_unannotated"].append(group.loc[(group['coltype'] == 'cell') & (group['chrscale'] == False) & (group['annotated'] == False), 'num_species'].values[0])
        yaxis_bot_species["total"].append(group.loc[(group['coltype'] == 'total'), 'num_species'].values[0])

    colorkey = {"chrscale_annotated": "green",
                "chrscale_unannotated": "blue",
                "nonchrscale_annotated": "orange",
                "nonchrscale_unannotated": "red",
                "total": "black"
                }

    # convert both of the xaxes to datetime
    for thisax in [ax1, ax2]:
        thisax.xaxis.set_major_locator(mdates.MonthLocator())
        thisax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))

    # make the plot, add a legend and axis labels
    # plot the genomes on the top panel
    for key in yaxis_top_genomes:
        ax1.plot(xaxis_plotdate, yaxis_top_genomes[key], label=key, color=colorkey[key])

    ax1.legend()
    ax1.set_ylabel('Number of Genomes')
    # plot the species on the bottom panel
    for key in yaxis_bot_species:
        ax2.plot(xaxis_plotdate, yaxis_bot_species[key], label=key, color=colorkey[key])
    ax2.legend()
    ax2.set_ylabel('Number of Species')

    # make the x-axis labels every 3 months from the earliest date to the latest date. Start in January of the earliest year.
    #  This is a bit hacky, but it works.
    #  The first date is the earliest date, the last date is the latest date
    # use the datetime package
    firstdate = xaxis_plotdate[0]
    lastdate  = xaxis_plotdate[-1]
    firstyear = firstdate.year
    lastyear  = lastdate.year+1
    # make a list of all the months from January to December, Just get every 3 months
    months = [str(x) for x in range(1, 13, 3)]
    years = [str(x) for x in range(firstyear, lastyear+1)]
    months_years = []
    for year in years:
        for month in months:
            months_years.append('{}-{}'.format(year, month))
    months_years_dt = [datetime.strptime(x, '%Y-%m') for x in months_years]
    months_years_str = [x.strftime('%Y-%m-%d') for x in months_years_dt]
    for thisax in [ax1, ax2]:
        thisax.set_xticks(months_years_dt)
        # set the xticklabels to be the months and years
        thisax.set_xticklabels(months_years_str, rotation=90, fontsize=6)

    # return the figure
    return fig

def main():
    args = parse_args()

    # load the input table as a tsv. The columns are coltype, chrscale, annotated, num_species, num_genomes, assembly_release_date_cutoff
    df = pd.read_csv(args.input, sep='\t')
    # parse the columns of chrscale and annotated as booleans, parsing the strings
    for col in ['chrscale', 'annotated']:
        df[col] = df[col].apply(lambda x: True if x == 'True' else False)
        df[col] = df[col].astype(bool)

    # if the path to the output pdf does not end in .pdf, add it
    if not args.output.endswith('.pdf'):
        args.output += '.pdf'
    # if the directory in which we want to save the pdf does not exist, make it safely, recursively
    outdir = os.path.dirname(args.output)
    if not os.path.isdir(outdir) and outdir != '':
        os.makedirs(outdir)
    # open a pdf with matplotlib for writing all the plots to
    with PdfPages(args.output) as pdf:
        # make the line plot of the number of genomes and species
        thisplt = lineplot_of_my_filtered_assemblies(df)
        # tight layout
        thisplt.tight_layout()
        pdf.savefig(thisplt)

    # print out the table
    print(df)

if __name__ == '__main__':
    main()