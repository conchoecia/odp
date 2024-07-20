#!/usr/bin/env python

"""
This performs fourier transforms on a set of rates from the plot_branch_stats_vs_time.py plot
"""

import argparse
import os

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from   scipy.signal import hann
from scipy.signal import correlate
from scipy.fftpack import fft, fftfreq
from scipy.signal import medfilt, correlate

from matplotlib.backends.backend_pdf import PdfPages

# we need to do this for plotting text
import odp_plotting_functions as odp_plot

def parse_args():
    """
    We need the path to a tsv file that contains the change rates over time.
      - rates:   the path to the tsv file containing the rates
      - agecol:  the name of the column containing the ages
      - ratecol: the name of the column containing the rates
      - minage: the minimum age to consider
      - polynomial: the degree of the polynomial to fit to the data
      - outprefix:  the prefix of the output files
      - custom_peaks: comma-separated list of peaks to plot
    """
    parser = argparse.ArgumentParser(description='Perform fourier transforms on a set of rates')
    parser.add_argument('--rates', type=str, help='Path to the tsv file containing the rates')
    parser.add_argument('--agecol', type=str, default='age', help='The name of the column containing the ages')
    parser.add_argument('--ratecol', type=str, default='fusion_rate_at_this_age_mean', help='The name of the column containing the rates')
    parser.add_argument('--minage', type=float, default=-1, help='The minimum age to consider')
    parser.add_argument('--polynomial', type=int, default=3, help='The degree of the polynomial to fit to the data')
    parser.add_argument('--outprefix', type=str, default='fourier', help='The prefix of the output files')
    parser.add_argument('--custom_peaks', type=str, default='', help='Comma-separated list of peaks (floats, periods) to plot')
    # check that the rates file exists
    args = parser.parse_args()
    if not os.path.exists(args.rates):
        parser.error(f'The file {args.rates} does not exist')
    if args.custom_peaks:
        args.custom_peaks = [float(peak) for peak in args.custom_peaks.replace(' ', '').split(',')]
    return args

def fft_unpadded(time, detrended_values):
    """
    Does a Fourier transform on the time series with x=time, y=values.
    """
    N = len(time)
    T = np.mean(np.diff(time))
    yf_original = fft(detrended_values)
    xf_original = fftfreq(N, T)[:N//2]
    spectrum = 2.0/N * np.abs(yf_original[:N//2])
    return spectrum, xf_original, yf_original

def fft_padded(time, values, zero_padding_factor=10):
    """
    Does a Fourier transform on the time series with x=time, y=values.
    """
    # Zero-pad the data
    N = len(time)
    T = np.mean(np.diff(time))
    padded_length = N * zero_padding_factor
    values_padded = np.pad(values, (0, padded_length - N), 'constant')

    # Perform Fourier transform on the zero-padded data
    yf_padded = fft(values_padded)
    xf_padded = fftfreq(padded_length, T)[:padded_length//2]
    padded_spectrum = 2.0/padded_length * np.abs(yf_padded[:padded_length//2])
    return padded_spectrum, xf_padded, yf_padded

def fourier_of_time_series(df, timecol, valuecol, polynomial, outprefix, zero_padding_factor=10, custom_peaks = []):
    """
    Perform a fourier transform on the time series with x=time, y=values.
    All of the x and y values will be ints or floats.

    The input variables:
      - df:  the dataframe containing the time series. The time is not necesasrily in order.
      - timecol: the column name of the time values. The values of the column should be evenly spaced.
      - valuecol: the column name of the values to be transformed.

    The dataframe will be sorted by the time column, then those two columns will be extracted for the
    Fourier transform. This will be plotted and saved to a dataframe.
    """
    outpdf = outprefix + '.pdf'
    # Drop NaN values and sort by time column
    df = df.dropna(subset=[valuecol, timecol]).sort_values(by=[timecol])
    # remove any row if time is 0
    df = df[df[timecol] != 0].reset_index(drop=True)
    time = df[timecol].values
    N = len(time)
    values = df[valuecol].values
    # ensure that the time values are evenly spaced
    time_diff = np.diff(time)
    print(df)
    if not np.allclose(time_diff, time_diff[0]):
        raise ValueError(f'The time values are not evenly spaced: {time_diff}')
    values = df[valuecol].values
    # fit a cubic function to the data
    p = np.polyfit(time, values, polynomial)
    detrended_values = values - np.polyval(p, time)

    # Perform Fourier transform on the original data
    mag_spectrum_original, xf_original, yf_original = fft_unpadded(time, detrended_values)

    # Zero-pad the data
    mag_spectrum_padded, xf_padded, yf_padded = fft_padded(time, detrended_values, zero_padding_factor=10)

    # PLOTTING
    # make three plots stacked on top of each other
    fig, axs = plt.subplots(3, 1, figsize=(10, 15))

    # Top panel: Original time series with cubic fit
    axs[0].plot(time, values, color='black', label='Original Data')
    axs[0].plot(time, np.polyval(p, time), 'b-', label=f'Cubic Fit (Degree {polynomial})')
    axs[0].set_xlabel('Time')
    axs[0].set_ylabel('Values')
    axs[0].legend()

    # Top middle panel: Fourier transform (zero-padded and original)
    axs[1].plot(xf_padded, mag_spectrum_padded, color='black', label='Zero-padded Spectrum')
    ax2 = axs[1].twinx()
    ax2.plot(xf_original, mag_spectrum_original, color='red', alpha=0.25)
    ax2.set_ylim(0, 1.1 * np.max(mag_spectrum_original))
    ax2.yaxis.label.set_color('red')
    ax2.tick_params(axis='y', colors='red')
    axs[1].plot([], [], color='red', label='Original Spectrum')
    axs[1].set_xlim(0, 0.1)
    axs[1].set_ylim(0, 1.1 * np.max(mag_spectrum_padded))
    axs[1].set_xlabel('Frequency')
    axs[1].set_ylabel('Magnitude (Zero-padded Spectrum)')
    ax2.set_ylabel('Magnitude (Original Spectrum)', color='red')

    default_peaks = np.argsort(np.abs(yf_original[:N//2]))[-3:]
    default_peaks = [xf_original[x] for x in default_peaks]
    # Bottom middl panel: Residuals and sine waves
    colorrange = ["#4F8FB8", "#FFA07A", "#90EE90"]
    if len(custom_peaks) > 0:
        peaks = custom_peaks
    else:
        peaks = default_peaks
    print(f"The peaks for {outpdf} are {peaks}")
    counter = 0
    for peak in peaks:
        period = 1/peak
        axs[1].axvline(peak, color=colorrange[counter], linestyle='--', label = f'Peak at {peak:.4f}, Period {period:.2f}')
        # just get the top of the box
        y_plotpos = 0.5 * np.max(mag_spectrum_padded)
        counter += 1
    axs[1].legend()

    # now we actually plot the new one
    axs[2].plot(time, detrended_values, color='black', label='Detrended Data')
    # Plot the top peaks as sine waves with optimal alignment
    for i, peak in enumerate(peaks):
        period = 1 / peak
        sine_wave = np.sin(2 * np.pi * time / period) * np.max(detrended_values) * 0.5

        # Compute optimal alignment using cross-correlation
        correlation = correlate(detrended_values, sine_wave, mode='full')
        lags = np.arange(-len(time) + 1, len(time))
        optimal_lag = lags[np.argmax(correlation)]
        shifted_sine_wave = np.sin(2 * np.pi * (time - time[optimal_lag]) / period) * np.max(detrended_values) * 0.5

        axs[2].plot(time, shifted_sine_wave, color=colorrange[i], linestyle='--', label=f'Sine Wave, Period {period:.2f}')

    axs[2].legend()
    axs[2].set_xlabel('Time')
    axs[2].set_ylabel('Detrended Values')

    # make the plot title be the outpdf
    axs[0].set_title(outpdf)
    # Save the plot to a PDF file
    plt.savefig(outpdf)
    plt.close()


    # Now run the rsimulations and make new PDFs
    rsims = r_simulations(time, detrended_values, default_peaks, mag_spectrum_original, polynomial, outprefix+"_unpadded", padded = False, simulations = 5000)
    rsims = r_simulations(time, detrended_values, peaks,         mag_spectrum_padded,   polynomial, outprefix+"_padded", padded = True, simulations = 5000)

def exponential_background(magnitudes, heights):
    """
    This function calculates the probability function P_h(f) for a given range of frequencies f.
    """
    h = magnitudes
    b = np.mean(heights)

    return 1 - np.exp(-h / b)

def r_simulations(time, values, peaks_to_test, magnitudes, polynomial, outprefix, padded = False, simulations = 30000):
    """
    This performs the "R" type Monte Carlo simulation from Rohde & Muller (2005) Nature.
    The description of this simulation, taken directly from the paper, is as follows:

    Choosing an appropriate model for a Monte Carlo simulation is always to some degree a matter of
      educated guesswork based on the expected behavior of the system. In our case, several authors
      have  proposed  that  the  changes  observed  in  diversity  are  effectively a random walk.
      So our “R” model for the changes in diversity is simply to  construct  random  walks  by  randomly
      rearranging  the steps  between  bins  in  the  existing  data.  This has the added advantage that
      diversity both starts and ends in the same place  in  every  simulation  and  will  generally  trend
      upward. After removing the best fitting cubic trend, the power  spectra  for  these  simulations  were
      computed  and  their  average  taken  to  obtain  the  R  background  given  in  the  paper. Random  walks
      heavily  favor  low  frequency changes, and this is reflected in the strong trend of the R background.
    """
    # get the average height over all frequencies
    magnitudes_average = np.mean(magnitudes)
    pdf_path = outprefix + '_rsims.pdf'
    # make a pdf for the r outs
    tsv_path = outprefix + '_rsims.tsv'
    tsv_entries = []
    with PdfPages(pdf_path) as pdf:
        chunk_to_support = []
        for num_chunks in range(10, 50):
            entries = []
            for i in range(simulations):
                #print(f"    Simulation {i}/{simulations}", end='\r')
                start = values[0]
                # get the transitions and randomize them
                #transitions = np.random.permutation(np.diff(values))
                # THIS IS TYPE W
                # calculate the transitions, split into chunks of length 27, then shuffle the chunks
                transitions = np.diff(values)
                chunks = np.array_split(transitions, num_chunks)
                # Shuffle the chunks
                shuffled_chunks = [chunks[z] for z in np.random.permutation(len(chunks))]
                # Flatten the shuffled chunks
                transitions = np.concatenate(shuffled_chunks)
                # for detrended_shuffled, randomize the transitions and add them to the start, then i=-1 sequentially
                shuffled = [start]
                for j in range(len(transitions)):
                    shuffled.append(shuffled[-1] + transitions[j])
                # ensure that the lengths are the same
                if len(values) != len(shuffled):
                    raise ValueError("The detrended_shuffled does not have the same length as the detrended_values")
                # enforce that the detrended_shuffled has the same start and stop as the detrended_values, rounded to 5 decimal places
                if not np.isclose(values[-1], shuffled[-1], atol=1e-5):
                    print(values[-1], shuffled[-1])
                    raise ValueError("The detrended_shuffled does not have the same start and stop as the detrended_values")
                # fit a polynomial to the shuffled data
                p = np.polyfit(time, shuffled, polynomial)
                detrended = shuffled - np.polyval(p, time)
                if padded:
                    # perform padded fourier transform
                    mag_spectrum_padded, xf_padded, yf_padded = fft_padded(time, detrended)
                else:
                    # perform unpadded fourier transform
                    mag_spectrum_padded, xf_padded, yf_padded = fft_unpadded(time, detrended)
                # get the mean magnitudes of the simulated peaks
                simulated_mag_means = np.mean(mag_spectrum_padded)
                # scale the mag_spectrum_padded s.t. the average magnitude is the same as the original
                mag_spectrum_padded = mag_spectrum_padded * magnitudes_average / simulated_mag_means
                # get the measurements at each of the peaks
                thisentry = {}
                for peak in peaks_to_test:
                    # get the value at this peak from mag_spectrum_padded and add it to thisentry
                    peak_index = np.argmin(np.abs(xf_padded - peak))
                    thisentry[peak] = mag_spectrum_padded[peak_index]
                entries.append(thisentry)
            # get the magnitudes of the real peaks
            real_magnitudes = {}
            for peak in peaks_to_test:
                peak_index = np.argmin(np.abs(xf_padded - peak))
                real_magnitudes[peak] = magnitudes[peak_index]
            valuesdf = pd.DataFrame(entries)
            print(valuesdf)

            # make a new page that has 1 row and n columns. n is the number of peaks
            fig, axs = plt.subplots(1, len(peaks_to_test), figsize=(5 * len(peaks_to_test), 5))
            # set the title to have the number of chunks and the number of simulations and the outprefix
            fig.suptitle(f"Chunks {num_chunks}, Simulations {simulations}, {outprefix}")
            peaks_support = {"chunks": num_chunks, "simulations": simulations, "outprefix": outprefix}
            for peak_i in range(len(peaks_to_test)):
                peak = peaks_to_test[peak_i]
                # for each peak, x is the max of the magnitude from the simulations
                x = np.linspace(0, np.max(valuesdf[peak]) , 1000)
                # get the number of values in the simulation that are less than x
                y = []
                for i in range(len(x)):
                    y.append(np.sum(valuesdf[peak] < x[i])/simulations)
                # convert the y-axis to log
                #y = np.log10(y)
                # create a dataframe and store it in dfs
                df = pd.DataFrame({'magnitude': x, 'percent': y})
                print("the peak is: ", peak)
                print("The period is: ", 1/peak)
                print("The chunk is: ", num_chunks)
                print("The magnitude of the peak is: ", real_magnitudes[peak])
                print("x% of the simulations are less than the real magnitude: ", np.sum(valuesdf[peak] < real_magnitudes[peak])/simulations)
                period = 1/peak
                peak_support = np.sum(valuesdf[peak] < real_magnitudes[peak])/simulations
                peaks_support[f"period_{period:.2f}_observed"] = peak_support
                # draw a horizontal, red dotted line, label with text
                axs[peak_i].axhline(peak_support, color='red', linestyle='--', label=f"Support: {peak_support:.3f}")
                # plot the percent on the y-axis, log scale, but labeled as a percentage
                # plot the magnitude on the x-axis
                # plot a vertical line at the observed axis
                # The color of the data will be blue, and the vertical line will be black
                axs[peak_i].plot(x, y, color='blue')
                axs[peak_i].axvline(real_magnitudes[peak], color='black')
                # set the title to be the period of the peak
                axs[peak_i].set_title(f"Period {1/peak:.2f}")
                # set the x-axis to be the magnitude
                axs[peak_i].set_xlabel('Magnitude')
                # set the y-axis to be the percent
                axs[peak_i].set_ylabel('Percent')
                # set the y-axis to be log of percentage
                # we want these percents to be plotted
                #y_ticks = [1, 60, 90, 96, 99, 99.6, 99.9, 99.96, 99.99, 99.996, 99.999]
                #y_ticks_log = np.log10([y/100 for y in y_ticks])
                #y_tick_labels = [f'{y:.3f}%' for y in y_ticks]
                #axs[peak_i].set_yticks(y_ticks_log)
                #axs[peak_i].set_yticklabels(y_tick_labels)

                # plot the probabilities
                probabilities_exp = exponential_background(df['magnitude'], valuesdf[peak])
                # find the magnitude index that is closest to the real magnitude
                real_magnitude_index = np.argmin(np.abs(df['magnitude'] - real_magnitudes[peak]))
                probabilities_exp_support = probabilities_exp[real_magnitude_index]
                peaks_support[f"period_{period:.2f}_expected"] = probabilities_exp_support
                # plot this as a green, finely dotted line
                axs[peak_i].plot(df['magnitude'], probabilities_exp, color='green', linestyle=':', label=f'Expected Support: {probabilities_exp_support:.3f}')
                # make a legend
                axs[peak_i].legend()
            chunk_to_support.append(peaks_support)
            pdf.savefig()
            plt.close()
        # make a plot of chunk_to_support
        chunk_to_support_df = pd.DataFrame(chunk_to_support)
        # plot the column "chunks" as x, and the columns beginning with "period_" as y.
        # These are n lines on the same figure. Include a legend.
        # there is a second panel that is a boxplot of the support columns
        fig, ax = plt.subplots(1, 2, figsize=(10, 5))
        fig.suptitle(f"Support for {outprefix}")
        support_cols = [col for col in chunk_to_support_df.columns if "period_" in col]
        for col in support_cols:
            ax[0].plot(chunk_to_support_df["chunks"], chunk_to_support_df[col], label=col)

        ax[0].legend()
        ax[0].set_xlabel("Chunks")
        ax[0].set_ylabel("Percent of Simulations\nLess Than Real Magnitude")

        changetype = ""
        if "loss" in outprefix:
            changetype = "losses"
        elif "fusions" in outprefix:
            changetype = "fusions"
        treatment = ""
        if "padded" in outprefix:
            treatment = "padded"
        elif "unpadded" in outprefix:
            treatment = "unpadded"
        # Clade	change_type	treatment	period	support_min	support_max	support_mean	support_median	support_std	polynomial
        # calculate the mean, median, and std of the support columns.
        # put the text on the figure
        outstring = ""
        for col in support_cols:
            obs_exp = ""
            if "observed" in col:
                obs_exp = "observed"
            elif "expected" in col:
                obs_exp = "expected"
            mean = chunk_to_support_df[col].mean()
            minv = chunk_to_support_df[col].min()
            maxv = chunk_to_support_df[col].max()
            median = chunk_to_support_df[col].median()
            std = chunk_to_support_df[col].std()
            outstring += f"{col} Min: {minv:.3f} Max: {maxv:.3f} Mean: {mean:.3f} Median: {median:.3f} Std: {std:.3f}\n"
            tsv_entries.append({"clade": outprefix, "obs_exp": obs_exp, "change_type": changetype, "treatment": treatment,
                                "period": float(col.split("_")[1]), "support_min": minv, "support_max": maxv, "support_mean": mean,
                                "support_median": median, "support_std": std, "polynomial": polynomial})
        ax[0].text(0+0.01, 1-0.01, outstring, fontsize=5, va = "top", transform=ax[0].transAxes)
        # set ylim to 0, 1
        ax[0].set_ylim(0, 1)
        ax[1].boxplot(chunk_to_support_df[support_cols])
        ax[1].set_ylim(0, 1)
        ax[1].set_xticklabels([col for col in support_cols], rotation=90, fontsize=5)
        ax[1].set_xlabel("Periods")
        ax[1].set_ylabel("Percent of Simulations\nLess Than Real Magnitude")
        # fit this to the pdf
        plt.tight_layout()
        pdf.savefig()
        plt.close()
        print(chunk_to_support_df)
        print("TSV entries is: ", tsv_entries)

    # convert the tsv_entries to a dataframe and save it to a tsv
    tsv_df = pd.DataFrame(tsv_entries)
    tsv_df.to_csv(tsv_path, sep='\t', index=False)

def main():
    args = parse_args()
    print(args)
    odp_plot.format_matplotlib()
    df = pd.read_csv(args.rates, sep='\t')
    # make sure that ratecol is in df.columns
    if args.ratecol not in df.columns:
        raise ValueError(f'The column {args.ratecol} is not in the dataframe')
    # make sure that agecol is in df.columns
    if args.agecol not in df.columns:
        raise ValueError(f'The column {args.agecol} is not in the dataframe')
    print(df.columns)
    if args.minage < 0:
        df = df[df[args.agecol] > args.minage]
    else:
        df = df[df[args.agecol] < args.minage]
    print("Args.custom_peaks", args.custom_peaks)
    if len(args.custom_peaks) > 0:
        fourier_of_time_series(df, args.agecol, args.ratecol, args.polynomial, args.outprefix, custom_peaks=args.custom_peaks)
    else:
        fourier_of_time_series(df, args.agecol, args.ratecol, args.polynomial, args.outprefix)

if __name__ == '__main__':
    main()