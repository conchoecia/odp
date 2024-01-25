#!/usr/bin/env python
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

#import os
#import sys
#snakefile_path = os.path.dirname(os.path.realpath(workflow.snakefile))
#dependencies_path = snakefile_path
#sys.path.insert(1, dependencies_path)
import perspchrom_df_to_tree as pdtt

configfile: "config.yaml"
config["sims_per_file"] = 50

rule all:
    input:
        expand("simresults/dfsim_run_{simnum}.tsv",
                simnum = range(int(config["numsims"]/config["sims_per_file"])))

rule sim:
    input:
        sampledf = config["sampledf"],
        rbhdf    = config["rbhdf"]
    output:
        simresultstsv = "simresults/dfsim_run_{simnum}.tsv"
    params:
        sims_per_file = config["sims_per_file"]
    threads: 1
    resources:
        runtime = 60,   # around 100 simulations takes 20 minutes
        mem_mb  = 1200
    run:
        pdtt.run_n_simulations_save_results(input.sampledf,
                                           input.rbhdf,
                                           output.simresultstsv,
                                           num_sims=params.sims_per_file,
                                           abs_bin_size=25,
                                           frac_bin_size=0.05,
                                           verbose = True)
