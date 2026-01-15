#!/usr/bin/env python

"""
This performs fourier transforms on a set of rates from the plot_branch_stats_vs_time.py plot
"""

import argparse
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.signal import correlate, find_peaks
from scipy.fftpack import fft, fftfreq
from scipy import stats

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
    parser.add_argument('--rates',        type=str, help='Path to the tsv file containing the rates')
    parser.add_argument('--agecol',       type=str, default='age', help='The name of the column containing the ages')
    parser.add_argument('--ratecol',      type=str, default='fusion_rate_at_this_age_mean', help='The name of the column containing the rates')
    parser.add_argument('--minage',       type=float, default=-1, help='[DEPRECATED] Use --min_time instead. The minimum age to consider')
    parser.add_argument('--min_time',     type=float, default=1.0, help='Minimum age in MYA for periodicity analysis (default: 1.0)')
    parser.add_argument('--max_time',     type=float, default=542.0, help='Maximum age in MYA for periodicity analysis (default: 542.0)')
    parser.add_argument('--polynomial',   type=int, default=3, help='The degree of the polynomial to fit to the data')
    parser.add_argument('--outprefix',    type=str, default='fourier', help='The prefix of the output files')
    parser.add_argument('--custom_peaks', type=str, default='', help='Comma-separated list of peaks (floats, periods) to plot')
    parser.add_argument('--include_unpadded', action='store_true', help='Also generate analysis for unpadded spectrum (default: padded only)')
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

def find_spectrum_peaks(spectrum, frequencies, n_peaks=3, use_peak_detection=True):
    """
    Find distinct peaks in the Fourier spectrum.
    
    Parameters:
    -----------
    spectrum : array
        The magnitude spectrum
    frequencies : array
        The corresponding frequencies
    n_peaks : int
        Number of peaks to return (default: 3)
    use_peak_detection : bool
        If True, use scipy.signal.find_peaks to find local maxima (good for smooth padded spectra).
        If False, simply take top N magnitudes (needed for monotonically decreasing unpadded spectra).
    
    Returns:
    --------
    peak_indices : array
        Indices of the peaks in the spectrum
    peak_frequencies : array
        Frequencies of the peaks (sorted by magnitude, descending)
    """
    if use_peak_detection:
        # Find all local maxima without any distance constraint
        peaks, properties = find_peaks(spectrum)
        
        # Sort peaks by magnitude (descending) and take top n_peaks
        if len(peaks) > 0:
            peak_magnitudes = spectrum[peaks]
            sorted_indices = np.argsort(peak_magnitudes)[::-1]  # descending order
            top_peaks = peaks[sorted_indices[:n_peaks]]
            top_frequencies = frequencies[top_peaks]
            return top_peaks, top_frequencies
        else:
            # Fallback: if no peaks found, use simple magnitude sorting
            use_peak_detection = False
    
    if not use_peak_detection:
        # Simple approach: take top N magnitudes (for monotonically decreasing spectra)
        top_indices = np.argsort(spectrum)[-n_peaks:][::-1]  # descending order
        top_frequencies = frequencies[top_indices]
        return top_indices, top_frequencies

def fourier_of_time_series(df, timecol, valuecol, polynomial, outprefix, zero_padding_factor=10, custom_peaks = [], include_unpadded=False):
    """
    Orchestrator function that runs Fourier analysis for padded spectrum (always) 
    and optionally for unpadded spectrum.
    
    Parameters:
    -----------
    df : DataFrame
        Input dataframe with time series data
    timecol : str
        Column name for time values
    valuecol : str
        Column name for values to transform
    polynomial : int
        Degree of polynomial for detrending
    outprefix : str
        Prefix for output files
    zero_padding_factor : int
        Factor for zero-padding (default: 10)
    custom_peaks : list
        Custom peak frequencies to use (default: [])
    include_unpadded : bool
        Whether to also generate unpadded analysis (default: False)
    """
    # Prepare data once
    df = df.dropna(subset=[valuecol, timecol]).sort_values(by=[timecol])
    df = df[df[timecol] != 0].reset_index(drop=True)
    time = df[timecol].values
    values = df[valuecol].values
    
    # Verify evenly spaced
    time_diff = np.diff(time)
    print(df)
    if not np.allclose(time_diff, time_diff[0]):
        raise ValueError(f'The time values are not evenly spaced: {time_diff}')
    
    # Always run padded analysis
    print("\n" + "="*60)
    print("RUNNING PADDED SPECTRUM ANALYSIS")
    print("="*60)
    run_single_analysis(time, values, polynomial, outprefix + "_padded", 
                       spectrum_type="padded", custom_peaks=custom_peaks)
    
    # Optionally run unpadded analysis
    if include_unpadded:
        print("\n" + "="*60)
        print("RUNNING UNPADDED SPECTRUM ANALYSIS")
        print("="*60)
        run_single_analysis(time, values, polynomial, outprefix + "_unpadded", 
                           spectrum_type="unpadded", custom_peaks=custom_peaks)

def run_single_analysis(time, values, polynomial, outprefix, spectrum_type, custom_peaks=[]):
    """
    Perform Fourier transform on the time series with x=time, y=values.
    All of the x and y values will be ints or floats.

    The input variables:
      - time: array of time values (evenly spaced)
      - values: array of values to be transformed
      - polynomial: degree of polynomial for detrending
      - outprefix: prefix for output files (should include _padded or _unpadded suffix)
      - spectrum_type: "padded" or "unpadded" - determines which spectrum to use for peak detection
      - custom_peaks: list of custom peak frequencies to use instead of auto-detection

    The function will create a PDF with the Fourier analysis results.
    """
    outpdf = outprefix + '.pdf'
    N = len(time)
    # fit a cubic function to the data
    p = np.polyfit(time, values, polynomial)
    detrended_values = values - np.polyval(p, time)

    # Perform Fourier transform on the original data
    mag_spectrum_original, xf_original, yf_original = fft_unpadded(time, detrended_values)

    # Zero-pad the data
    mag_spectrum_padded, xf_padded, yf_padded = fft_padded(time, detrended_values, zero_padding_factor=10)

    # Detect peaks from the appropriate spectrum based on spectrum_type
    if spectrum_type == "padded":
        # Use find_peaks for smooth padded spectrum
        default_peak_indices, default_peaks_freq = find_spectrum_peaks(
            mag_spectrum_padded, xf_padded, n_peaks=3, use_peak_detection=True
        )
        print(f"Detected peaks from PADDED spectrum at frequencies: {default_peaks_freq}")
    else:  # unpadded
        # Use simple magnitude sorting for potentially monotonic unpadded spectrum
        default_peak_indices, default_peaks_freq = find_spectrum_peaks(
            mag_spectrum_original, xf_original, n_peaks=3, use_peak_detection=False
        )
        print(f"Detected peaks from UNPADDED spectrum at frequencies: {default_peaks_freq}")
    
    print(f"Corresponding periods (MYA): {[1/f for f in default_peaks_freq]}")
    
    if len(custom_peaks) > 0:
        peaks_freq = custom_peaks
    else:
        peaks_freq = default_peaks_freq
    
    # Run simulations before creating figure to determine correct panel count
    support_threshold = 0.85
    min_chunks = 20
    n_significant = 0
    
    cache_check = load_cached_simulation_data(outprefix)
    if cache_check is not None:
        print("Found cached simulation data - using for layout calculation")
        chunk_to_support = cache_check
        # Calculate n_significant from cache
        for peak in peaks_freq:
            period_val = 1/peak
            support_col = f"period_{period_val:.2f}_observed"
            if support_col in cache_check.columns:
                chunk_mask = cache_check["chunks"] > min_chunks
                filtered_data = cache_check[chunk_mask][support_col]
                if len(filtered_data) > 0 and filtered_data.median() > support_threshold:
                    n_significant += 1
        print(f"Pre-computed {n_significant} significant periods for figure layout")
    else:
        # No cache - run simulations now to determine layout before creating figure
        print(f"No cache found - running {spectrum_type} simulations to determine layout...")
        if spectrum_type == "padded":
            chunk_to_support = r_simulations(time, detrended_values, peaks_freq, mag_spectrum_padded, 
                                            polynomial, outprefix, padded=True, simulations=5000)
        else:  # unpadded
            chunk_to_support = r_simulations(time, detrended_values, peaks_freq, mag_spectrum_original, 
                                            polynomial, outprefix, padded=False, simulations=5000)
        
        # Calculate n_significant from fresh simulation results
        for peak in peaks_freq:
            period_val = 1/peak
            support_col = f"period_{period_val:.2f}_observed"
            if support_col in chunk_to_support.columns:
                chunk_mask = chunk_to_support["chunks"] > min_chunks
                filtered_data = chunk_to_support[chunk_mask][support_col]
                if len(filtered_data) > 0 and filtered_data.median() > support_threshold:
                    n_significant += 1
        print(f"Calculated {n_significant} significant periods from new simulations")
    
    # PLOTTING
    # Create figure with base 4 panels + validated comparison (if significant) + sequential removal panels
    # If n_significant > 0: add 1 panel for validated oscillatory fit comparison
    total_panels = 4 + (1 if n_significant > 0 else 0) + n_significant
    fig, axs = plt.subplots(total_panels, 1, figsize=(12, 5 * total_panels))

    # Panel A: Original time series with polynomial fit AND detrended residuals
    axs[0].plot(time, values, color='black', linewidth=1, label='Original Data')
    axs[0].plot(time, np.polyval(p, time), 'b-', linewidth=1, label=f'Polynomial Fit (Degree {polynomial})')
    axs[0].plot(time, detrended_values, 'r-', linewidth=1, alpha=0.7, label='Detrended Residuals')
    axs[0].set_xlabel('Time (MYA)')
    axs[0].set_ylabel('Rate / Residuals')
    axs[0].legend()
    axs[0].set_title('Panel A: Raw Data + Polynomial Fit + Residuals')

    # Use peaks detected earlier
    colorrange = ["#4F8FB8", "#FFA07A", "#90EE90"]
    peaks = peaks_freq
    print(f"The peaks for {outpdf} are {peaks}")

    # Panel B: Fourier spectrum (zero-padded and original)
    axs[1].plot(xf_padded, mag_spectrum_padded, color='black', linewidth=1, label='Zero-padded Spectrum')
    ax2 = axs[1].twinx()
    ax2.plot(xf_original, mag_spectrum_original, color='red', linewidth=1, alpha=0.25)
    ax2.set_ylim(0, 1.1 * np.max(mag_spectrum_original))
    ax2.yaxis.label.set_color('red')
    ax2.tick_params(axis='y', colors='red')
    axs[1].plot([], [], color='red', label='Original Spectrum')
    axs[1].set_xlim(0, 0.1)
    axs[1].set_ylim(0, 1.1 * np.max(mag_spectrum_padded))
    axs[1].set_xlabel('Frequency (cycles per MYR)')
    axs[1].set_ylabel('Magnitude (Zero-padded Spectrum)')
    ax2.set_ylabel('Magnitude (Original Spectrum)', color='red')
    
    # Title indicates which spectrum was used for peak detection
    spectrum_source = "Padded" if spectrum_type == "padded" else "Unpadded"
    axs[1].set_title(f'Panel B: Fourier Spectrum (Peaks from {spectrum_source})')
    
    # Add period labels at the top inside the plot
    # Create custom period tick locations at meaningful frequencies
    period_values = [500, 250, 140, 100, 62, 50, 40, 30, 25, 20]  # periods in MYA
    freq_values = [1.0/p for p in period_values if 1.0/p <= 0.1]  # convert to frequencies within xlim
    period_labels = [str(int(1.0/f)) for f in freq_values]
    
    # Add text annotations for periods at the top of the plot
    y_pos = 0.95  # Position near top of plot in axis coordinates
    for freq, label in zip(freq_values, period_labels):
        axs[1].text(freq, y_pos, label, transform=axs[1].get_xaxis_transform(),
                    ha='center', va='top', fontsize=8, rotation=0)
    
    # Add a label for the period axis
    axs[1].text(0.5, 0.98, 'Period (MYA)', transform=axs[1].transAxes,
                ha='center', va='top', fontsize=9, weight='bold')
    
    # Mark peaks on spectrum
    counter = 0
    for peak in peaks:
        period = 1/peak
        axs[1].axvline(peak, color=colorrange[counter], linestyle='--', 
                      label=f'Peak at {peak:.4f}, Period {period:.2f} MYA')
        counter += 1
    
    # Move legend to lower right to avoid overlap with period labels
    axs[1].legend(fontsize=8, loc='lower right')

    # Panel C: Residuals with linear and oscillatory fits
    # Fit linear regression to residuals
    slope, intercept, r_value_linear, p_value_linear, std_err = stats.linregress(time, detrended_values)
    linear_fit = slope * time + intercept
    r_squared_linear = r_value_linear ** 2
    
    # Debug output for linear regression
    print(f"Linear regression debug:")
    print(f"  slope={slope:.6e}, intercept={intercept:.6e}")
    print(f"  R²={r_squared_linear:.6f}, p={p_value_linear:.6e}")
    print(f"  Time range: {time.min():.2f} to {time.max():.2f}")
    print(f"  Residuals: mean={np.mean(detrended_values):.6e}, std={np.std(detrended_values):.6e}")
    print(f"  Residuals range: {detrended_values.min():.6e} to {detrended_values.max():.6e}")
    
    # Plot residuals
    axs[2].plot(time, detrended_values, color='black', linewidth=1, alpha=0.5, label='Residuals')
    axs[2].axhline(0, color='gray', linestyle=':', linewidth=0.5)
    
    # Plot linear fit in red
    axs[2].plot(time, linear_fit, 'r-', linewidth=2, label=f'Linear Regression Fit (R²={r_squared_linear:.4f}, p={p_value_linear:.2e})')
    
    # Build composite oscillatory fit from all peaks
    
    # Build composite oscillatory fit from all peaks
    oscillatory_fit = np.zeros_like(time, dtype=float)
    for i, peak in enumerate(peaks):
        period = 1 / peak
        sine_wave = np.sin(2 * np.pi * time / period) * np.max(detrended_values) * 0.5
        # Compute optimal alignment using cross-correlation
        correlation = correlate(detrended_values, sine_wave, mode='full')
        lags = np.arange(-len(time) + 1, len(time))
        optimal_lag = lags[np.argmax(correlation)]
        shifted_sine_wave = np.sin(2 * np.pi * (time - time[optimal_lag]) / period) * np.max(detrended_values) * 0.5
        oscillatory_fit += shifted_sine_wave
    
    # Normalize oscillatory fit to match residual scale
    oscillatory_fit = oscillatory_fit / len(peaks)
    
    # Compute R² for oscillatory fit
    ss_res_osc = np.sum((detrended_values - oscillatory_fit) ** 2)
    ss_tot = np.sum((detrended_values - np.mean(detrended_values)) ** 2)
    r_squared_osc = 1 - (ss_res_osc / ss_tot)
    
    # Build multi-line legend label with periods
    periods_str = ', '.join([f'{1/p:.1f}' for p in peaks])
    osc_label = f'Composite Sine Wave Fit (R²={r_squared_osc:.4f})\nPeriods: {periods_str} MYA'
    
    # Plot oscillatory fit in blue
    axs[2].plot(time, oscillatory_fit, 'b-', linewidth=2, label=osc_label)
    
    axs[2].set_xlabel('Time (MYA)')
    axs[2].set_ylabel('Residuals')
    axs[2].set_title('Panel C: Residuals with Linear vs Oscillatory Fits')
    axs[2].legend(fontsize=8)
    
    # Add text box with comparison
    if r_squared_osc > r_squared_linear:
        comparison_text = f'Oscillatory fit is better\n(ΔR² = {r_squared_osc - r_squared_linear:.4f})'
    else:
        comparison_text = f'Linear fit is better\n(ΔR² = {r_squared_linear - r_squared_osc:.4f})'
    axs[2].text(0.02, 0.98, comparison_text, transform=axs[2].transAxes, 
                fontsize=9, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    # Panel D: Simulation summary - box-and-whisker plots
    # (chunk_to_support already populated from cache or fresh simulations above)
    support_cols_obs = [col for col in chunk_to_support.columns if "period_" in col and "observed" in col]
    support_cols_exp = [col for col in chunk_to_support.columns if "period_" in col and "expected" in col]
    
    print(f"\n{spectrum_type.capitalize()} box plot columns:")
    print(f"  Found {len(support_cols_obs)} observed columns: {support_cols_obs}")
    print(f"  Found {len(support_cols_exp)} expected columns: {support_cols_exp}")
    print(f"  All columns: {list(chunk_to_support.columns)}")
    
    # Create pairs (Expected, Observed) sorted by period (ascending = shortest first)
    periods = sorted(set([float(col.split('_')[1]) for col in support_cols_obs]))
    data_to_plot = []
    labels = []
    for period in periods:
        # Add Expected first
        exp_col = f"period_{period:.2f}_expected"
        if exp_col in support_cols_exp:
            data_to_plot.append(chunk_to_support[exp_col].values)
            labels.append(f'{period:.1f}\nExp')
        # Then Observed
        obs_col = f"period_{period:.2f}_observed"
        if obs_col in support_cols_obs:
            data_to_plot.append(chunk_to_support[obs_col].values)
            labels.append(f'{period:.1f}\nObs')
    
    bp = axs[3].boxplot(data_to_plot, labels=labels, patch_artist=True)
    # Color Expected boxes green, Observed boxes blue (alternating pattern)
    for i, box in enumerate(bp['boxes']):
        if i % 2 == 0:  # Even indices are Expected
            box.set_facecolor('lightgreen')
        else:  # Odd indices are Observed
            box.set_facecolor('lightblue')
    
    axs[3].set_ylabel("Support (% Sims < Real)", fontsize=10)
    axs[3].set_ylim(0, 1)
    axs[3].set_title(f"Panel D: {spectrum_type.capitalize()} Simulation Support (Observed vs Expected)", fontsize=11)
    axs[3].tick_params(axis='x', labelsize=8)
    axs[3].axhline(0.95, color='red', linestyle='--', alpha=0.3, linewidth=1)
    axs[3].text(0.98, 0.95, '95%', transform=axs[3].get_yaxis_transform(), 
                ha='right', va='bottom', fontsize=8, color='red')

    # Add sequential period removal panels to main figure (if significant periods exist)
    if n_significant > 0:
        # Filter for significant periods using simulation results
        significant_periods = []
        for peak in peaks_freq:
            period_val = 1/peak
            support_col = f"period_{period_val:.2f}_observed"
            if support_col in chunk_to_support.columns:
                chunk_mask = chunk_to_support["chunks"] > min_chunks
                filtered_data = chunk_to_support[chunk_mask][support_col]
                if len(filtered_data) > 0:
                    median_support = filtered_data.median()
                    if median_support > support_threshold:
                        significant_periods.append((period_val, median_support, peak))
        
        # Sort by support descending
        significant_periods.sort(key=lambda x: x[1], reverse=True)
        
        print(f"\nAdding validated oscillatory fit comparison and {len(significant_periods)} sequential removal panels:")
        for period_val, support, freq in significant_periods:
            print(f"  Period {period_val:.1f} MYA: median support = {support:.3f}")
        
        # Panel E: Linear vs Oscillatory fit using ONLY significant peaks
        axs[4].plot(time, detrended_values, color='black', linewidth=1, alpha=0.5, label='Residuals')
        axs[4].axhline(0, color='gray', linestyle=':', linewidth=0.5)
        
        # Plot linear fit
        axs[4].plot(time, linear_fit, 'r-', linewidth=2, label=f'Linear Fit (R²={r_squared_linear:.4f})')
        
        # Build composite oscillatory fit from ONLY significant peaks
        validated_oscillatory_fit = np.zeros_like(time, dtype=float)
        for i, (period_val, support, freq) in enumerate(significant_periods):
            sine_wave = np.sin(2 * np.pi * time / period_val) * np.max(detrended_values) * 0.5
            # Compute optimal alignment using cross-correlation
            correlation = correlate(detrended_values, sine_wave, mode='full')
            lags = np.arange(-len(time) + 1, len(time))
            optimal_lag = lags[np.argmax(correlation)]
            shifted_sine_wave = np.sin(2 * np.pi * (time - time[optimal_lag]) / period_val) * np.max(detrended_values) * 0.5
            validated_oscillatory_fit += shifted_sine_wave
        
        # Normalize oscillatory fit to match residual scale
        validated_oscillatory_fit = validated_oscillatory_fit / len(significant_periods)
        
        # Compute R² for validated oscillatory fit
        ss_res_val_osc = np.sum((detrended_values - validated_oscillatory_fit) ** 2)
        ss_tot = np.sum((detrended_values - np.mean(detrended_values)) ** 2)
        r_squared_val_osc = 1 - (ss_res_val_osc / ss_tot)
        
        # Build legend label with validated periods
        validated_period_labels = [f"{p:.1f} MYA" for p, s, f in significant_periods]
        validated_legend_label = f'Validated Oscillatory Fit (R²={r_squared_val_osc:.4f})\nPeriods: ' + ', '.join(validated_period_labels)
        
        axs[4].plot(time, validated_oscillatory_fit, 'b-', linewidth=2, label=validated_legend_label)
        axs[4].set_xlabel('Time (MYA)')
        axs[4].set_ylabel('Residuals')
        axs[4].set_title('Panel E: Residuals with Linear vs Validated Oscillatory Fits')
        axs[4].legend(fontsize=8)
        
        # Add comparison text box
        if r_squared_val_osc > r_squared_linear:
            comparison_text = f'Validated oscillatory fit is better\n(ΔR² = {r_squared_val_osc - r_squared_linear:.4f})'
        else:
            comparison_text = f'Linear fit is better\n(ΔR² = {r_squared_linear - r_squared_val_osc:.4f})'
        axs[4].text(0.02, 0.98, comparison_text, transform=axs[4].transAxes, 
                    fontsize=9, verticalalignment='top',
                    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        current_residuals = detrended_values.copy()
        colorrange_seq = ["#4F8FB8", "#FFA07A", "#90EE90", "#FFD700", "#FF69B4"]
        
        for i, (period_val, support, freq) in enumerate(significant_periods):
            panel_idx = 5 + i
            panel_letter = chr(70 + i)  # F, G, H, etc.
            
            # Plot current residuals
            axs[panel_idx].plot(time, current_residuals, 'k-', linewidth=1, label='Residuals')
            axs[panel_idx].axhline(0, color='gray', linestyle=':', linewidth=0.5)
            
            # Fit sine wave to current residuals
            sine_wave = fit_sine_wave_to_residuals(time, current_residuals, period_val)
            
            color = colorrange_seq[i % len(colorrange_seq)]
            axs[panel_idx].plot(time, sine_wave, color=color, linestyle='--', linewidth=2, 
                           label=f'Period {period_val:.1f} MYA (Support: {support:.3f})')
            
            axs[panel_idx].set_xlabel('Time (MYA)', fontsize=10)
            axs[panel_idx].set_ylabel('Residuals', fontsize=10)
            axs[panel_idx].set_title(f'Panel {panel_letter}: Sequential Removal Step {i+1} - Removing {period_val:.1f} MYA Period', fontsize=11)
            axs[panel_idx].legend(fontsize=9)
            axs[panel_idx].grid(alpha=0.3)
            
            # Subtract this sine wave for next iteration
            current_residuals = current_residuals - sine_wave
    else:
        print("\nNo significant periods found - no sequential removal panels added")

    # Overall title and save
    fig.suptitle(outpdf, fontsize=14, y=0.998)
    plt.tight_layout(rect=[0, 0, 1, 0.995])
    plt.savefig(outpdf)
    plt.close()

def exponential_background(magnitudes, heights):
    """
    This function calculates the probability function P_h(f) for a given range of frequencies f.
    """
    h = magnitudes
    b = np.mean(heights)

    return 1 - np.exp(-h / b)

def fit_sine_wave_to_residuals(time, residuals, period):
    """
    Fit a sine wave with optimal phase alignment to residuals.
    Returns the fitted sine wave.
    """
    # Create sine wave
    sine_wave = np.sin(2 * np.pi * time / period)
    
    # Compute optimal alignment using cross-correlation
    correlation = correlate(residuals, sine_wave, mode='full')
    lags = np.arange(-len(time) + 1, len(time))
    optimal_lag = lags[np.argmax(correlation)]
    
    # Shift sine wave for optimal alignment
    if optimal_lag >= 0 and optimal_lag < len(time):
        shifted_sine_wave = np.sin(2 * np.pi * (time - time[optimal_lag]) / period)
    else:
        shifted_sine_wave = sine_wave
    
    # Scale amplitude to match residuals using least squares
    # Find amplitude that minimizes sum of squared residuals
    amplitude = np.sum(residuals * shifted_sine_wave) / np.sum(shifted_sine_wave ** 2)
    
    return amplitude * shifted_sine_wave

def load_cached_simulation_data(outprefix):
    """
    Check if cached simulation data exists and load it.
    Returns the chunk_to_support dataframe if found, None otherwise.
    """
    chunk_support_path = outprefix + '_chunk_support.tsv'
    if os.path.exists(chunk_support_path):
        print(f"Found cached simulation data: {chunk_support_path}")
        return pd.read_csv(chunk_support_path, sep='\t')
    return None

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
        if "_padded" in outprefix:
            treatment = "padded"
        elif "_unpadded" in outprefix:
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
    
    # Save chunk_to_support_df for caching/reuse
    chunk_support_path = outprefix + '_chunk_support.tsv'
    chunk_to_support_df.to_csv(chunk_support_path, sep='\t', index=False)
    print(f"Saved chunk support data to: {chunk_support_path}")
    
    # Return the summary dataframe for use in main plot
    return chunk_to_support_df

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
    
    # Show age range in data before filtering
    age_min = df[args.agecol].min()
    age_max = df[args.agecol].max()
    print(f"Age range in data: {age_min:.2f} - {age_max:.2f} (n={len(df)} rows)")
    
    # Handle deprecated --minage argument
    if args.minage >= 0:
        print(f"WARNING: --minage is deprecated. Use --min_time and --max_time instead.")
        print(f"Filtering to ages >= {args.minage} Mya (ignoring --min_time and --max_time)")
        df = df[df[args.agecol] >= args.minage]
    else:
        # Apply min_time and max_time filters (default: 1-542 MYA)
        # Ages are encoded as negative (0=present, -1=1MYA, -542=542MYA)
        # So filter to ages between -max_time and -min_time
        print(f"Filtering to time range: {args.min_time} - {args.max_time} MYA")
        print(f"  (filtering to age values: {-args.max_time} to {-args.min_time})")
        df = df[(df[args.agecol] >= -args.max_time) & (df[args.agecol] <= -args.min_time)]
    
    print(f"After filtering: {len(df)} rows remaining")
    
    # Check if dataframe is empty after filtering
    if len(df) == 0:
        print(f"\nERROR: No data remaining after filtering to {args.min_time}-{args.max_time} MYA")
        print(f"Original data age range was {age_min:.2f} - {age_max:.2f}")
        print(f"Adjust --min_time and --max_time to match your data range.")
        sys.exit(1)
    
    print("Args.custom_peaks", args.custom_peaks)
    if len(args.custom_peaks) > 0:
        fourier_of_time_series(df, args.agecol, args.ratecol, args.polynomial, args.outprefix, 
                              custom_peaks=args.custom_peaks, include_unpadded=args.include_unpadded)
    else:
        fourier_of_time_series(df, args.agecol, args.ratecol, args.polynomial, args.outprefix,
                              include_unpadded=args.include_unpadded)

if __name__ == '__main__':
    main()