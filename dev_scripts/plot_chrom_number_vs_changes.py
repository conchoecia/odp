#!/usr/bin/env python

"""
This script plots the number of chromosomes in a sample vs the number of changes.
The user must provide to
sys.arvg[1] the path to the 'per_species_ALG_presence_fusions.tsv' file
sys.arvg[2] the path to the 'species_chrom_counts.tsv' file
"""

# odp stuff to format the plot
import os
import sys
# ODP-specific imports
thisfile_path = os.path.dirname(os.path.realpath(__file__))
scripts_path = os.path.join(thisfile_path, "../scripts")
sys.path.insert(1, scripts_path)
import odp_plotting_functions as odp_plot

# matplotlib stuff
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, TwoSlopeNorm

from perspchrom_df_to_tree import parse_gain_loss_from_perspchrom_df
import pandas as pd
import numpy as np
import scipy
from scipy.special import comb
from scipy.stats import gaussian_kde
import sys

def panel_chromnum_vs_fusions_hexbin(ax, labels, chromnum, colocs, losses):
    """
    Hexbin plot of chromosome number vs fusion count to address overplotting
    """
    hb = ax.hexbin(chromnum, colocs, gridsize=40, cmap='Blues', 
                   bins='log', mincnt=1, linewidths=0.2, edgecolors='lightgray')
    ax.set_xlim(0, 100)
    ax.set_xlabel("Number of chromosomes")
    ax.set_ylabel("Number of colocalizations (fusions)")
    ax.set_title("Chromosome count vs Fusions (hexbin density)")
    
    # Add colorbar
    cbar = plt.colorbar(hb, ax=ax)
    cbar.set_label('log10(count)', rotation=270, labelpad=15)
    
    # Add correlation stats
    spearman_corr, p_value = scipy.stats.spearmanr(chromnum, colocs)
    kendall_tau, kp_value = scipy.stats.kendalltau(chromnum, colocs)
    ax.text(0.02, 0.98, "Spearman: {:.3f} (p={:.2e})\nKendall: {:.3f} (p={:.2e})".format(
        spearman_corr, p_value, kendall_tau, kp_value), 
        transform=ax.transAxes, va='top', fontsize=9,
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    return ax

def panel_chromnum_vs_losses_hexbin(ax, labels, chromnum, colocs, losses):
    """
    Hexbin plot of chromosome number vs loss count (absolute values)
    """
    losses_abs = [-1 * x for x in losses]  # Convert to positive
    hb = ax.hexbin(chromnum, losses_abs, gridsize=40, cmap='Reds', 
                   bins='log', mincnt=1, linewidths=0.2, edgecolors='lightgray')
    ax.set_xlim(0, 100)
    ax.set_xlabel("Number of chromosomes")
    ax.set_ylabel("Number of losses (absolute)")
    ax.set_title("Chromosome count vs Losses (hexbin density)")
    
    # Add colorbar
    cbar = plt.colorbar(hb, ax=ax)
    cbar.set_label('log10(count)', rotation=270, labelpad=15)
    
    # Add correlation stats
    spearman_corr, p_value = scipy.stats.spearmanr(chromnum, losses_abs)
    kendall_tau, kp_value = scipy.stats.kendalltau(chromnum, losses_abs)
    ax.text(0.02, 0.98, "Spearman: {:.3f} (p={:.2e})\nKendall: {:.3f} (p={:.2e})".format(
        spearman_corr, p_value, kendall_tau, kp_value), 
        transform=ax.transAxes, va='top', fontsize=9,
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    return ax

def panel_coloc_vs_losses_hexbin(ax, labels, chromnum, colocs, losses):
    """
    Hexbin plot of colocalizations vs losses to show correlation density
    """
    losses_abs = [-1 * x for x in losses]  # Convert to positive
    hb = ax.hexbin(colocs, losses_abs, gridsize=35, cmap='Greys', 
                   bins='log', mincnt=1, linewidths=0.2, edgecolors='lightgray')
    ax.set_xlabel("Number of colocalizations")
    ax.set_ylabel("Number of losses (absolute)")
    ax.set_title("Fusions vs Losses (hexbin density)")
    
    # Add colorbar
    cbar = plt.colorbar(hb, ax=ax)
    cbar.set_label('log10(count)', rotation=270, labelpad=15)
    
    # Add correlation stats
    spearman_corr, p_value = scipy.stats.spearmanr(colocs, losses_abs)
    kendall_tau, kp_value = scipy.stats.kendalltau(colocs, losses_abs)
    ax.text(0.5, 0.2, "Spearman: {:.3f} (p={:.2e})\nKendall: {:.3f} (p={:.2e})".format(
        spearman_corr, p_value, kendall_tau, kp_value), 
        transform=ax.transAxes, fontsize=9,
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    return ax

def panel_fraction_fusions_and_loss_counts(ax, labels, chromnum, colocs, losses, num_ALGs=29):
    """
    Stacked panel: fraction of possible fusions (top) and absolute loss counts (bottom)
    """
    # Calculate fraction of possible fusions
    fraction_coloc = []
    chromnum_frac = []
    for i in range(len(chromnum)):
        num_combinations = comb(num_ALGs + losses[i], 2)
        if num_combinations > 0:
            fraction_coloc.append(colocs[i] / num_combinations)
            chromnum_frac.append(chromnum[i])
    
    # Get absolute loss values
    losses_abs = [-1 * x for x in losses]
    
    # Plot fraction of fusions in top half (positive y)
    ax.scatter(chromnum_frac, fraction_coloc, color='blue', 
               alpha=0.15, s=20, edgecolors='none', label='Fraction fusions')
    
    # Plot absolute loss counts in bottom half (negative y for visual separation)
    losses_neg = [-1 * x for x in losses_abs]  # Make negative for bottom half
    ax.scatter(chromnum, losses_neg, color='red', 
               alpha=0.15, s=20, edgecolors='none', label='Loss counts')
    
    # Add horizontal line at y=0
    ax.axhline(0, color='black', linewidth=1.5, alpha=0.7)
    
    # Add reference line at 29 chromosomes
    ymin = min(losses_neg) if losses_neg else -30
    ymax = max(fraction_coloc) if fraction_coloc else 0.5
    ax.plot([29, 29], [ymin, ymax], color='gray', linewidth=1, 
            alpha=0.5, linestyle='--')
    ax.text(29, ymax * 0.95, '29 BCnS ALGs', rotation=90, 
            va='top', ha='right', fontsize=9, alpha=0.7)
    
    ax.set_xlim(0, 100)
    ax.set_xlabel('Number of chromosomes')
    ax.set_ylabel('Fraction fusions (top) | Loss counts (bottom, negative)')
    ax.set_title('Fusion fraction vs Loss counts')
    ax.legend(loc='upper right', fontsize=9)
    
    # Add correlation stats for each half
    if len(chromnum_frac) > 0:
        sp_frac, pv_frac = scipy.stats.spearmanr(chromnum_frac, fraction_coloc)
        ax.text(0.02, 0.98, 'Fusions: ρ={:.3f}'.format(sp_frac),
                transform=ax.transAxes, va='top', fontsize=8,
                bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.7))
    
    sp_loss, pv_loss = scipy.stats.spearmanr(chromnum, losses_abs)
    ax.text(0.02, 0.02, 'Losses: ρ={:.3f}'.format(sp_loss),
            transform=ax.transAxes, va='bottom', fontsize=8,
            bbox=dict(boxstyle='round', facecolor='lightcoral', alpha=0.7))
    
    return ax

def panel_chromsize_vs_changes_bins(ax, labels, chromnum, colocs, losses):
    """
    Put the changes into bins instead of doing a scatter plot
    """
    fusion_bins = {} # keys are tuples: (number of chromosomes, number of fusions)
    loss_bins   = {} # keys are tuples: (number of chromosomes, number of losses)

    for i in range(len(chromnum)):
        # colocalizations
        fusion_bin = (chromnum[i], colocs[i])
        if fusion_bin not in fusion_bins:
            fusion_bins[fusion_bin] = 0
        fusion_bins[fusion_bin] += 1
        # losses
        loss_bin   = (chromnum[i], losses[i])
        if loss_bin not in loss_bins:
            loss_bins[loss_bin] = 0
        loss_bins[loss_bin] += 1

    #print("fusion bins are: {}".format(fusion_bins))
    max_divisor = 10
    # make a color scale of the number of samples in each bin. 0 is white, max is blue. Max is he max value of colocs
    colormap = plt.cm.Blues
    # get the max number of colocs when the second value of the bin is not 0
    max_fusions_nonzero = max([fusion_bins[x] for x in fusion_bins if x[1] != 0])/max_divisor
    norm = Normalize(vmin=0, vmax=max_fusions_nonzero)

    # now make rectangles for each coloc bin
    for thisbin in fusion_bins:
        # get the color for this bin
        color = colormap(norm(fusion_bins[thisbin]))
        # get the x and y coordinates of the bin
        x = thisbin[0]
        y = thisbin[1]
        # make a rectangle for this bin
        rect = mpatches.Rectangle((x,y), 1, 1,
                                  linewidth=0, edgecolor='none',
                                  facecolor=color)
        ax.add_patch(rect)

    colormap2 = plt.cm.Reds
    max_losses_nonzero = max([loss_bins[x] for x in loss_bins if x[1] != 0])/max_divisor
    norm2 = Normalize(vmin=0, vmax=max_losses_nonzero)

    # now make rectangles for each coloc bin
    for thisbin in loss_bins:
        # get the color for this bin
        color = colormap2(norm2(loss_bins[thisbin]))
        # get the x and y coordinates of the bin
        x = thisbin[0]
        y = thisbin[1]
        # make a rectangle for this bin
        rect = mpatches.Rectangle((x,y-1), 1, 1,
                                  linewidth=0, edgecolor='none',
                                  facecolor=color)
        ax.add_patch(rect)

    max_dimension = 100
    # change the x-axis to be between 0 and 100
    ax.set_xlim(0, max_dimension)
    # change the y-axis to be between 0 and the max number of fusions
    ax.set_ylim(-30, -30 + max_dimension)

    # make a thin vertical line at 29 - the number of chromosomes in the human genome
    ax.plot([29, 29], [-30, 100], color="black", lw=1, alpha=0.5)
    # label this vertically as "29 BCnS ALGs"
    ax.text(29, -20, "29 BCnS ALGs", rotation=90, verticalalignment="center", horizontalalignment="right")

    # x-axis label
    ax.set_xlabel("Number of chromosomes")
    ax.set_ylabel("number of changes. colocalizations are blue, losses are negative\n Color scale maxes out at 1/10 the max count")

    return ax

def panel_chromsize_vs_changes(ax, labels, chromnum, colocs, losses):
    """
    generates a single panel for the plot of the number of chromosomes vs the number of changes
    """

    # On the panel plot the number of fusions (y1) on the positive y axis
    #  and the number of losses (y2) on the negative y axis. They both use x.
    ax.scatter(chromnum, colocs, color="blue", lw = 0, alpha = 0.02)
    ax.scatter(chromnum, losses, color="red",  lw = 0, alpha = 0.02)

    # make a black vertical line between the two points
    for i in range(len(chromnum)):
        ax.plot([chromnum[i], chromnum[i]],
                [colocs[i], losses[i]], color="black", lw=2, alpha=0.01)

    # add a legend of the two types of dots
    # make the dots' alpha 1 so they are visible
    ax.scatter([], [], color="blue", label="fusions", lw = 0, alpha = 1)
    ax.scatter([], [], color="red", label="losses", lw = 0, alpha = 1)
    ax.legend(loc="upper right")

    # make the x-axis go between 0 and 100
    ax.set_xlim(0, 100)
    # make the x-axis label
    ax.set_xlabel("Number of chromosomes")
    ax.set_ylabel("number of changes. colocalizations are blue, losses are negative")
    ax.set_title("Number of chromosomes vs number of changes")
    return ax

def panel_chromosomes_vs_losses(ax, labels, chromnum, colocs, losses):
    """
    plots the number of chromosomes(x) vs the number of losses(y)
    """
    ax.scatter(chromnum, losses, color="black", lw = 0, alpha = 0.05)
    ax.set_xlim(0, np.median(chromnum)*3)
    ax.set_xlabel("Number of chromosomes")
    ax.set_ylabel("Number of losses")

    # Calculate Spearman's rank correlation coefficient
    spearman_corr, p_value = scipy.stats.spearmanr(chromnum, losses)

    # Calculate Kendall's tau
    kendall_tau, kp_value = scipy.stats.kendalltau(chromnum, losses)

    ax.text(0.5, 0.2, "Spearman's coeff.: {:.3f}\nP-value: {:e}\n\nKendall's tau coeff: {:.3f}\nP-value: {:e}".format(
        spearman_corr, p_value, kendall_tau, kp_value), transform=ax.transAxes)

    print(f"chromnum vs. losses Kendall's tau: {kendall_tau}")
    print(f"chromnum vs. losses P-value: {p_value}")

    print(f"chromnum vs. losses Spearman's rank correlation coefficient: {spearman_corr}")
    print(f"chromnum vs. losses P-value: {kp_value}")

    return ax

def panel_colocalization_vs_losses(ax, labels, chromnum, colocs, losses):
    """
    This panel is an x-y of the number of colocalizations vs the number of losses
    """
    # make a pandas df of colocs and losses
    df = pd.DataFrame({"colocs": colocs, "losses": losses})

    # On the panel plot the number of fusions (y1) on the positive y axis
    #  and the number of losses (y2) on the negative y axis. They both use x.
    ax.scatter(df["colocs"], df["losses"], color="black", lw = 0, alpha = 0.05)

    colocs = df["colocs"]
    losses = df["losses"]

    # Calculate Spearman's rank correlation coefficient
    spearman_corr, p_value = scipy.stats.spearmanr(colocs, losses)

    # Calculate Kendall's tau
    kendall_tau, kp_value = scipy.stats.kendalltau(colocs, losses)

    print(f"Kendall's tau: {kendall_tau}")
    print(f"P-value: {p_value}")

    print(f"Spearman's rank correlation coefficient: {spearman_corr}")
    print(f"P-value: {kp_value}")
    # Put the spearman's rank correlation coefficient on the plot. round to 3 decimal places
    # Can you put it in the lower right.
    ax.text(0.5, 0.2, "Spearman's coeff.: {:.3f}\nP-value: {:e}\n\nKendall's tau coeff: {:.3f}\nP-value: {:e}".format(
        spearman_corr, p_value, kendall_tau, kp_value), transform=ax.transAxes)

    # add axis labels
    ax.set_xlabel("Number of colocalizations")
    ax.set_ylabel("Number of losses")
    # make the dots' alpha 1 so they are visible
    ax.set_title("Correlation between number of colocalizations and losses")
    return ax

def panel_rank_coloc_losses(ax, x, y, xlabel, ylabel):
    """
    Plots the ranks of the number of colocalizations vs the number of losses
    """
    # Plotting ranks
    rank_x = scipy.stats.rankdata(x)
    rank_y = scipy.stats.rankdata(y)

    ax.scatter(rank_x, rank_y)
    ax.set_xlabel("Rank of {}".format(xlabel))
    ax.set_ylabel("Rank of {}".format(ylabel))
    ax.set_title("Rank of number of {} vs rank of {}".format(xlabel, ylabel))
    return ax

def panel_chromosomes_vs_fractionloss(ax, labels, chromnum, colocs, losses, num_ALGs = 29):
    """
    Makes a plot testing the spearman correlation of the number of chromosomes vs the fraction of losses.
    Don't count the ones with no losses
    """
    fraction_coloc = []
    new_chromnum   = []

    ignored_fraction_coloc = []
    ignored_chromnum       = []

    for i in range(len(chromnum)):
        # FIRST DO COLOC
        num_combinations = comb(num_ALGs + losses[i], 2)
        if num_combinations == 0:
            pass
            #fraction_coloc.append(0)
        else:
            thisfrac = colocs[i]/num_combinations
            if thisfrac == 0:
                ignored_fraction_coloc.append(thisfrac)
                ignored_chromnum.append(chromnum[i])
            else:
                fraction_coloc.append(thisfrac)
                new_chromnum.append(chromnum[i])

    ax.scatter(new_chromnum, fraction_coloc,
               color = "black", lw = 0, alpha = 0.05)
    ax.scatter(ignored_chromnum, ignored_fraction_coloc,
               color = "red", lw = 0, alpha = 0.05)
    ax.set_xlim(0, np.median(chromnum)*3)
    ax.set_xlabel("Number of chromosomes")
    ax.set_ylabel("Fraction of possible fusions")

    # Calculate Spearman's rank correlation coefficient
    spearman_corr, p_value = scipy.stats.spearmanr(new_chromnum, fraction_coloc)
    # Calculate Kendall's tau
    kendall_tau, kp_value = scipy.stats.kendalltau(new_chromnum, fraction_coloc)
    ax.text(0.5, 0.97, "All points\nSpearman's coeff.: {:.3f}\nP-value: {:e}\nKendall's tau coeff: {:.3f}\nP-value: {:e}".format(
        spearman_corr, p_value, kendall_tau, kp_value), va = "top", transform=ax.transAxes)

    # What is the correlation for samples with fewer than 20 chromosomes
    fraction_coloc_ltet20 = []
    new_chromnum_ltet20   = []
    for i in range(len(new_chromnum)):
        if new_chromnum[i] <= 20:
            fraction_coloc_ltet20.append(fraction_coloc[i])
            new_chromnum_ltet20.append(new_chromnum[i])

    thisx = new_chromnum_ltet20
    thisy = fraction_coloc_ltet20
    # Calculate Spearman's rank correlation coefficient
    spearman_corr, p_value = scipy.stats.spearmanr(thisx, thisy)
    # Calculate Kendall's tau
    kendall_tau, kp_value = scipy.stats.kendalltau(thisx, thisy)
    ax.text(0.5, 0.75, "chromnum <= 20\nSpearman's coeff.: {:.3f}\nP-value: {:e}\nKendall's tau coeff: {:.3f}\nP-value: {:e}".format(
        spearman_corr, p_value, kendall_tau, kp_value), va= "top", transform=ax.transAxes)

    # Make a legend
    ax.scatter([], [], label = "coeff. measured",
               color = "black", lw = 0, alpha = 0.8)
    ax.scatter([], [], label = "coeff. ignored",
               color="red", lw = 0, alpha = 0.8)
    ax.legend(loc="center right")

    return ax

def panel_fraction_fusions_losses_possible(ax, labels, chromnum, colocs, losses, num_ALGs = 29):
    """
    Makes a plot where the number of fusions and losses are plotted as a fraction of the potential.
    For example, if all 29 ALGs are present, there are 29 choose 2 possible combinations. We take the
    number of fusions and divide by the number of possible fusions. We do the same for losses, but # losses/29
    """
    # you can change this value if you want
    binsize = 0.0125

    loss_binsize = 1/num_ALGs
    coloc_bins = {} # keys are tuples: (number of chromosomes, fraction of fusions)
    loss_bins  = {}

    for i in range(len(chromnum)):
        # FIRST DO COLOC
        num_combinations = comb(num_ALGs + losses[i], 2)
        # print("There are {} colocs, {} losses, and {} possible combinations".format(colocs[i], losses[i], num_combinations))
        thisbin = (-1,-1)
        if num_combinations == 0:
            thisbin = (chromnum[i], float(0))
        else:
            frac_colocs = colocs[i]/num_combinations
            frac_bin = int(frac_colocs/binsize)*binsize
            if frac_bin == 0:
                frac_bin = float(0)
            thisbin = (chromnum[i], frac_bin)
        if thisbin not in coloc_bins:
            coloc_bins[thisbin] = 0
        coloc_bins[thisbin] += 1

        # NEXT DO LOSSES
        frac_losses = losses[i]/num_ALGs
        frac_bin = int(frac_losses/loss_binsize)*binsize # yes, we use binsize to scale to the top half
        thisbin = (chromnum[i], frac_bin)
        if thisbin not in loss_bins:
            loss_bins[thisbin] = 0
        loss_bins[thisbin] += 1

    max_divisor = 10
    # make a color scale of the number of samples in each bin. 0 is white, max is blue. Max is he max value of colocs
    colormap = plt.cm.Blues
    # get the max number of colocs when the second value of the bin is not 0
    max_fusions_nonzero = max([coloc_bins[x] for x in coloc_bins if x[1] != 0])/max_divisor
    norm = Normalize(vmin=0, vmax=max_fusions_nonzero)

    # now make rectangles for each coloc bin
    for thisbin in coloc_bins:
        # get the color for this bin
        color = colormap(norm(coloc_bins[thisbin]))
        # get the x and y coordinates of the bin
        x = thisbin[0]
        y = thisbin[1]
        # make a rectangle for this bin
        rect = mpatches.Rectangle((x,y), 1, binsize,
                                  linewidth=0, edgecolor='none',
                                  facecolor=color)
        ax.add_patch(rect)

    # now make rectangles for the loss bins
    colormap2 = plt.cm.Reds
    max_losses_nonzero = max([loss_bins[x] for x in loss_bins if x[1] != 0])/max_divisor
    norm2 = Normalize(vmin=0, vmax=max_losses_nonzero)
    # now make rectangles for each coloc bin
    for thisbin in loss_bins:
        # get the color for this bin
        color = colormap2(norm2(loss_bins[thisbin]))
        # get the x and y coordinates of the bin
        x = thisbin[0]
        y = thisbin[1]
        # make a rectangle for this bin
        rect = mpatches.Rectangle((x,y-binsize), 1, binsize,
                                  linewidth=0, edgecolor='none',
                                  facecolor=color)
        ax.add_patch(rect)

    ax.set_xlim(0, 100)
    ymin = min([x[1] for x in loss_bins])-binsize

    ymax = 0.5
    ax.set_ylim(ymin, ymax)
    # set the y-axis ticks to start at 0 and go to ymax
    ax.set_yticks(np.arange(0, ymax, 0.05))

    ax.set_ylabel("Fraction of possible fusions (top) and losses (bottom)")
    ax.set_xlabel("Number of chromosomes")

    return ax

def plot_chrom_number_vs_changes(changesfilename, chromsizefilename, outfilename):
    """
    saves a pdf of the plot of the number of chromosomes vs the number of changes
    """
    # CALL THIS TO GET THE VISUAL STYLE WE NEED
    odp_plot.format_matplotlib()

    # make a dictionary of the sample id to the number of chromosomes
    df = pd.read_csv(chromsizefilename, sep='\t')
    chromsize_to_chromnumber = dict(zip(df["sample"], df["chromosomes"]))
    del df

    changedf = parse_gain_loss_from_perspchrom_df(pd.read_csv(changesfilename, sep='\t'))
    samples_in_changedf = set(changedf["samplename"])

    labels   = []
    chromnum = [] # number of chromosomes
    colocs   = [] # number of fusions
    losses   = [] # number of losses
    for sample in chromsize_to_chromnumber:
        if (sample in samples_in_changedf) and (sample in chromsize_to_chromnumber):
            labels.append(sample)
            chromnum.append(chromsize_to_chromnumber[sample])
            subdf = changedf[changedf["samplename"] == sample]
            # y1 is the total size of the lists in the colocalizations column
            colocs.append(sum([len(x) for x in subdf["colocalizations"]]))
            # y2 is the total size of the lists in the losses column
            y2val = sum([len(x) for x in subdf["losses"]])
            losses.append(0 if y2val == 0 else -1*y2val)
            # if the number of losses is 0, print the sample name
            if y2val == 0:
                print(sample, " has 0 losses")

    fw = 26  # Increased width for 4th column
    fh = 38  # Increased height for additional row
    fig = plt.figure(figsize=(fw, fh))
    axes = []

    #for aligning all the panels
    left1   = 0.6
    left2   = 6.5
    left3   = 12.5
    left4   = 18.5

    axes = []

    bottom1 = 0.6
    paneldim = 5

    # This panel is the number of chromosomes vs the number of changes
    # we start with a single panel
    plot_params = [left1   /fw, # left offset
                   bottom1 /fh, # bottom offset
                   paneldim/fw, # width
                   paneldim/fh] # height
    axes.append(fig.add_axes(plot_params))
    axes[-1] = panel_chromsize_vs_changes(axes[-1], labels, chromnum, colocs, losses)

    # this panel is x-y of number of colocalizations vs number of losses
    plot_params = [left2   /fw, # left offset
                   bottom1 /fh, # bottom offset
                   paneldim/fw, # width
                   paneldim/fh] # height
    axes.append(fig.add_axes(plot_params))
    axes[-1] = panel_colocalization_vs_losses(axes[-1], labels, chromnum, colocs, losses)

    # add a panel of the rank plot of the numver of colocalizations and the number of losses
    plot_params = [left3   /fw, # left offset
                   bottom1 /fh, # bottom offset
                   paneldim/fw, # width
                   paneldim/fh] # height
    axes.append(fig.add_axes(plot_params))
    axes[-1] = panel_rank_coloc_losses(axes[-1], colocs, losses, "# of colocalizations", "# of losses")

    bottom2 = 6.5
    ## add a panel where we bin the number of chromosomes in each cell
    plot_params = [left1   /fw, # left offset
                   bottom2 /fh, # bottom offset
                   paneldim/fw, # width
                   paneldim/fh] # height
    axes.append(fig.add_axes(plot_params))
    axes[-1] = panel_chromsize_vs_changes_bins(axes[-1], labels, chromnum, colocs, losses)

    ## add a panel showing the number of chromsomes vs the number of losses
    plot_params = [left2   /fw, # left offset
                   bottom2 /fh, # bottom offset
                   paneldim/fw, # width
                   paneldim/fh] # height
    axes.append(fig.add_axes(plot_params))
    axes[-1] = panel_chromosomes_vs_losses(axes[-1], labels, chromnum, colocs, losses)

    # add a panel with the ranks of the number of chromosomes and the number of losses
    plot_params = [left3   /fw, # left offset
                   bottom2 /fh, # bottom offset
                   paneldim/fw, # width
                   paneldim/fh] # height
    axes.append(fig.add_axes(plot_params))
    axes[-1] = panel_rank_coloc_losses(axes[-1], chromnum, losses, "# of chromosomes", "# of losses")

    bottom3 = 12.25
    # add a panel where we bin the number of chromosomes in each cell
    plot_params = [left1   /fw, # left offset
                   bottom3 /fh, # bottom offset
                   paneldim/fw, # width
                   paneldim/fh] # height
    axes.append(fig.add_axes(plot_params))
    axes[-1] = panel_fraction_fusions_losses_possible(axes[-1], labels, chromnum, colocs, losses)

    # add a panel where we show the spearman correlation of fraction of fusions vs num_chromosomes
    plot_params = [left2   /fw, # left offset
                   bottom3 /fh, # bottom offset
                   paneldim/fw, # width
                   paneldim/fh] # height
    axes.append(fig.add_axes(plot_params))
    axes[-1] = panel_chromosomes_vs_fractionloss(axes[-1], labels, chromnum, colocs, losses)

    # NEW PANELS - Using 4th column and additional row for density visualizations
    
    # Row 1, Column 4: Hexbin for chromosome vs fusions
    plot_params = [left4   /fw, # left offset
                   bottom1 /fh, # bottom offset
                   paneldim/fw, # width
                   paneldim/fh] # height
    axes.append(fig.add_axes(plot_params))
    axes[-1] = panel_chromnum_vs_fusions_hexbin(axes[-1], labels, chromnum, colocs, losses)

    # Row 2, Column 4: Hexbin for chromosome vs losses
    plot_params = [left4   /fw, # left offset
                   bottom2 /fh, # bottom offset
                   paneldim/fw, # width
                   paneldim/fh] # height
    axes.append(fig.add_axes(plot_params))
    axes[-1] = panel_chromnum_vs_losses_hexbin(axes[-1], labels, chromnum, colocs, losses)

    # Row 3, Column 4: Hexbin for colocalizations vs losses
    plot_params = [left4   /fw, # left offset
                   bottom3 /fh, # bottom offset
                   paneldim/fw, # width
                   paneldim/fh] # height
    axes.append(fig.add_axes(plot_params))
    axes[-1] = panel_coloc_vs_losses_hexbin(axes[-1], labels, chromnum, colocs, losses)

    # Row 4: New row for the stacked fraction fusions (top) and loss counts (bottom)
    bottom4 = 18.0
    plot_params = [left1   /fw, # left offset
                   bottom4 /fh, # bottom offset
                   paneldim/fw, # width
                   paneldim/fh] # height
    axes.append(fig.add_axes(plot_params))
    axes[-1] = panel_fraction_fusions_and_loss_counts(axes[-1], labels, chromnum, colocs, losses)

    # Row 4, Column 2: Another view - 2D histogram version of fraction fusions
    plot_params = [left2   /fw, # left offset
                   bottom4 /fh, # bottom offset
                   paneldim/fw, # width
                   paneldim/fh] # height
    axes.append(fig.add_axes(plot_params))
    # Create 2D histogram of fraction fusions
    fraction_coloc_all = []
    chromnum_all = []
    for i in range(len(chromnum)):
        num_combinations = comb(29 + losses[i], 2)
        if num_combinations > 0:
            fraction_coloc_all.append(colocs[i] / num_combinations)
            chromnum_all.append(chromnum[i])
    
    if len(chromnum_all) > 0:
        h = axes[-1].hist2d(chromnum_all, fraction_coloc_all, bins=[50, 40], 
                           cmap='Blues', cmin=1)
        axes[-1].set_xlabel('Number of chromosomes')
        axes[-1].set_ylabel('Fraction of possible fusions')
        axes[-1].set_title('Fusion fraction density (2D histogram)')
        axes[-1].set_xlim(0, 100)
        cbar = plt.colorbar(h[3], ax=axes[-1])
        cbar.set_label('Count', rotation=270, labelpad=15)

    fig.savefig(outfilename, bbox_inches="tight")

def main():
    plot_chrom_number_vs_changes(sys.argv[1], sys.argv[2], "chrom_number_vs_changes.pdf")

if __name__ == "__main__":
    main()