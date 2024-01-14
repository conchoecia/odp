#!/usr/bin/env python

import perspchrom_df_to_tree as cdf

"""
Program  : perspchrom_df_stats_and_mc.snakefile
Language : snakemake
Date     : 2024-01-13
Author   : Darrin T. Schultz
Email    : darrin.schultz@univie.ac.at
Github   : https://github.com/conchoecia/odp
License  : GNU GENERAL PUBLIC LICENSE, Version 3, 29 June 2007. See the LICENSE file.

This snakefile controls the Monte Carlo analysis of testing the observed/expected
 ratio of chromosome fusion sizes. The python file that contains the code to perform
 this analysis is in the file perspchrom_df_to_tree.py. The reason that the MC is controlled
 with this script is the easy parallelization offered by snakemake and the integration
 with SLURM clusters.
"""

rule all:
    input:
        expand("simresults/dfsim_run_{simnum}.tsv",
                simnum = range(int(config["numsims"]/100)))

rule simulations:
    input:
        sampledf = config["sampledf"],
        rbhdf    = config["rbhdf"]
    output:
        simresultstsv = "simresults/dfsim_run_{simnum}.tsv"
    threads: 1
    resources:
        runtime = 20   # around 100 simulations takes 20 minutes
        mem_mb  = 1000
    run:
        cdf.run_n_simulations_save_results(input.sampledf,
                                           input.rbhdf,
                                           output.simresultstsv,
                                           num_sims=100,
                                           abs_bin_size=25,
                                           frac_bin_size=0.05)