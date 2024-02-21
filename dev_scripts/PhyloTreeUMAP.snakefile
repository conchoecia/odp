"""
This is a snakefile that is used to parallelize constructing a chromosomal linkage UMAP.
 For this UMAP, each dot is a genome or an inferred pseudogenome at one node.

I am writing a snakefile to do this, because it takes too long in a for loop.
"""

from PhyloTreeUMAP import (algcomboix_file_to_dict,
                           construct_coo_matrix_from_sampledf,
                           rbh_to_samplename,
                           rbh_directory_to_distance_matrix,
                           rbh_to_distance_gbgz,
                           ALGrbh_to_algcomboix,
                           plot_umap_from_files,
                           sampleToRbhFileDict_to_sample_matrix)
import os
import pandas as pd
from scipy.sparse import coo_matrix, lil_matrix, save_npz, load_npz

# Ignore all ResourceWarnings
import warnings
warnings.filterwarnings("ignore", category=ResourceWarning)

configfile: "config.yaml"
# First we define all of the RBH files
# get the rbh files in the directory
config["ALGname"] = "BCnSSimakov2022"

# CHECKS
# check that rbh_directory exists in the config file
if not "rbh_directory" in config:
    raise ValueError("rbh_directory not in config file")
# make sure that the user has provided the rbh file
if not "ALG_rbh_file" in config:
    raise ValueError("ALG_rbh_file not in config file")

# check that the ALG rbh file exists
if not os.path.exists(config["rbh_directory"]):
    raise ValueError(f"rbh_directory {config['rbh_directory']} does not exist")

config["rbh_directory"] = os.path.abspath(config["rbh_directory"])

config["rbh_files"] = list(sorted([os.path.join(config["rbh_directory"], f)
                           for f in os.listdir(config["rbh_directory"])
                           if f.endswith('.rbh')], reverse = True))
config["sample_to_rbh_file"] = {rbh_to_samplename(x, config["ALGname"]): x
                                for x in config["rbh_files"]}

# Results_base_directory is the prefix to which everything will be saved
results_base_directory = "GTUMAP"
if results_base_directory.endswith("/"):
    results_base_directory = results_base_directory[:-1]

rule all:
    input:
        expand(results_base_directory + "/AllSamples/AllSamples.neighbors_{n}.mind_{m}.missing_{sizeNaN}.df",
                n = [5, 10, 20, 50, 100, 250],
                m = [0.0, 0.1, 0.2, 0.5, 0.75, 0.9, 1.0],
                sizeNaN = ["small", "large"]),
        expand(results_base_directory + "/AllSamples/AllSamples.neighbors_{n}.mind_{m}.missing_{sizeNaN}.bokeh.html",
                sample = config["sample_to_rbh_file"].keys(),
                n = [5, 10, 20, 50, 100, 250],
                m = [0.0, 0.1, 0.2, 0.5, 0.75, 0.9, 1.0],
                sizeNaN = ["small", "large"]),
        #expand(results_base_directory + "/distance_matrices/{sample}.gb.gz",
        #                      sample = config["sample_to_rbh_file"].keys()),
        #results_base_directory + "/combo_to_index.txt",
        #results_base_directory + "/sampledf.tsv",


def generic_get_mem_mb(wildcards, attempt):
    """
    The amount of RAM needed for the script depends on the size of the input genome.
    """
    attemptdict = {1: 200,
                   2: 400,
                   3: 800,
                   4: 1600,
                   5: 3200,
                   6: 6400,
                  }
    return attemptdict[attempt]

rule samples_and_gzipped:
    input:
        rbh_file = lambda wildcards: config["sample_to_rbh_file"][wildcards.sample]
    output:
        gbgz = results_base_directory + "/distance_matrices/{sample}.gb.gz"
    threads: 1
    params:
        ALGname = config["ALGname"]
    retries: 6
    resources:
        mem_mb = generic_get_mem_mb,
        runtime = 2
    run:
        rbh_to_distance_gbgz(input.rbh_file, output.gbgz, params.ALGname)

rule sample_df:
    input:
        rbh_files = [config["sample_to_rbh_file"][x]
                     for x in config["sample_to_rbh_file"]]
    output:
        sampletsv = results_base_directory + "/sampledf.tsv"
    retries: 6
    threads: 1
    params:
        gbgz_directory = results_base_directory + "/distance_matrices"
    resources:
        mem_mb = 1000,
        runtime = 10
    run:
        sampleToRbhFileDict_to_sample_matrix(config["sample_to_rbh_file"],
                                             config["ALGname"],
                                             params.gbgz_directory,
                                             output.sampletsv)

rule combo_to_index:
    input:
        rbhfile = config["ALG_rbh_file"]
    output:
        outfile = results_base_directory + "/combo_to_index.txt"
    threads: 1
    resources:
        mem_mb = 1000,
        runtime = 5
    run:
        # We will need to calculate the rbh combo to index dictionary.
        # DO NOT bother reading in the existing file. It takes 5x longer
        #  to read in the file than it does to generate it.
        alg_combo_to_ix = ALGrbh_to_algcomboix(input.rbhfile)
        # save the dictionary pairs
        with open(output.outfile, "w") as f:
            for key, value in alg_combo_to_ix.items():
                f.write(f"{key}\t{value}\n")

def coo_get_mem_mb(wildcards, attempt):
    """
    The amount of RAM needed for the script depends on the size of the input genome.
    """
    attemptdict = {1: 10000,
                   2: 50000,
                   3: 100000,
                   4: 200000,
                   5: 300000,
                  }
    return attemptdict[attempt]

rule generate_coo_matrix:
    """
    This takes up a lot of RAM and time. For 3600 genomes, it took:
        CPU Efficiency: 97.32% of 01:57:32 core-walltime
        Job Wall-clock time: 01:57:32
        Memory Utilized: 216.61 GB
        Memory Efficiency: 110.91% of 195.31 GB
    """
    input:
        gbgzfiles = expand(results_base_directory + "/distance_matrices/{sample}.gb.gz",
                           sample = config["sample_to_rbh_file"].keys()),
        combotoindex = results_base_directory + "/combo_to_index.txt",
        sampletsv    = results_base_directory + "/sampledf.tsv"
    output:
        coo    = results_base_directory + "/allsamples.coo.npz"
    threads: 1
    retries: 5
    resources:
        mem_mb = coo_get_mem_mb,
        runtime = 300
    run:
        # read in the combo_to_index file as a df, then convert to a dict
        alg_combo_to_ix = algcomboix_file_to_dict(input.combotoindex)
        sampledf = pd.read_csv(input.sampletsv, sep = "\t", index_col = 0)
        coo = construct_coo_matrix_from_sampledf(sampledf, alg_combo_to_ix)
        save_npz(output.coo, coo)

rule plot_umap_of_files:
    input:
        sampletsv    = results_base_directory + "/sampledf.tsv",
        combotoindex = results_base_directory + "/combo_to_index.txt",
        coo    = results_base_directory + "/allsamples.coo.npz"
    output:
        df   = results_base_directory + "/AllSamples/AllSamples.neighbors_{n}.mind_{m}.missing_{sizeNaN}.df",
        html = results_base_directory + "/AllSamples/AllSamples.neighbors_{n}.mind_{m}.missing_{sizeNaN}.bokeh.html"
    threads: 1
    params:
        outdir = results_base_directory + "/all_species",
    retries: 5
    resources:
        mem_mb = coo_get_mem_mb,
        runtime = 60
    run:
        plot_umap_from_files(input.sampletsv, input.combotoindex, input.coo,
                             params.outdir, wildcards.sample, wildcards.sizeNaN, wildcards.n, wildcards.m)
