#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft, fftfreq
from scipy.signal import medfilt

# Example data (replace with your actual data)
time = np.linspace(0, 549, 550)
values = np.sin(2 * np.pi * time / 100) + 0.5 * np.random.normal(size=550)

# Perform the Fourier Transform
N = len(time)
T = 1  # Sample spacing (1 million years)
yf = fft(values)
xf = fftfreq(N, T)[:N//2]

# Compute the magnitude spectrum
mag_spectrum = 2.0/N * np.abs(yf[:N//2])

# Smooth the magnitude spectrum (e.g., median filter)
smoothed_spectrum = medfilt(mag_spectrum, kernel_size=5)

# Estimate the spectral background (e.g., using a rolling mean)
background_estimate = np.mean(smoothed_spectrum)

# Subtract the background from the original spectrum
clean_spectrum = mag_spectrum - background_estimate

# Plotting
plt.figure(figsize=(10, 6))

plt.subplot(2, 1, 1)
plt.plot(xf, mag_spectrum, label='Original Spectrum')
plt.plot(xf, smoothed_spectrum, label='Smoothed Spectrum')
plt.xlabel('Frequency (1/Millions of Years)')
plt.ylabel('Magnitude')
plt.title('Original and Smoothed Spectra')
plt.legend()
plt.grid(True)

plt.subplot(2, 1, 2)
plt.plot(xf, clean_spectrum, label='Clean Spectrum')
plt.xlabel('Frequency (1/Millions of Years)')
plt.ylabel('Magnitude')
plt.title('Spectral Background Subtraction')
plt.legend()
plt.grid(True)

plt.tight_layout()
# save this as a pdf
plt.savefig('spectral_background.pdf')
