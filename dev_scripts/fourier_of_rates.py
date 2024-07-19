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

def fourier_of_time_series(df, timecol, valuecol, polynomial, outpdf, zero_padding_factor=10, custom_peaks = []):
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
    # Drop NaN values and sort by time column
    df = df.dropna(subset=[valuecol, timecol]).sort_values(by=[timecol])
    # remove any row if time is 0
    df = df[df[timecol] != 0].reset_index(drop=True)
    time = df[timecol].values
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
    N = len(time)
    T = np.mean(np.diff(time))
    yf_original = fft(detrended_values)
    xf_original = fftfreq(N, T)[:N//2]
    mag_spectrum_original = 2.0/N * np.abs(yf_original[:N//2])

    # Zero-pad the data
    padded_length = N * zero_padding_factor
    detrended_values_padded = np.pad(detrended_values, (0, padded_length - N), 'constant')

    # Perform Fourier transform on the zero-padded data
    yf_padded = fft(detrended_values_padded)
    xf_padded = fftfreq(padded_length, T)[:padded_length//2]
    mag_spectrum_padded = 2.0/padded_length * np.abs(yf_padded[:padded_length//2])

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

    # Bottom middl panel: Residuals and sine waves
    colorrange = ["#4F8FB8", "#FFA07A", "#90EE90"]
    if len(custom_peaks) > 0:
        peaks = custom_peaks
    else:
        peaks = np.argsort(np.abs(yf_original[:N//2]))[-3:]
        peaks = [xf_original[x] for x in peaks]
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
    outpdf = args.outprefix + '.pdf'
    print("Args.custom_peaks", args.custom_peaks)
    if len(args.custom_peaks) > 0:
        fourier_of_time_series(df, args.agecol, args.ratecol, args.polynomial, outpdf, custom_peaks=args.custom_peaks)
    else:
        fourier_of_time_series(df, args.agecol, args.ratecol, args.polynomial, outpdf)

if __name__ == '__main__':
    main()