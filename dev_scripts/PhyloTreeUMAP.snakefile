"""
This is a snakefile that is used to parallelize constructing a chromosomal linkage UMAP.
 For this UMAP, each dot is a genome or an inferred pseudogenome at one node.

I am writing a snakefile to do this, because it takes too long in a for loop.
"""

from PhyloTreeUMAP import (algcomboix_file_to_dict,
                           construct_coo_matrix_from_sampledf,
                           filter_sample_df_by_clades,
                           rbh_to_samplename,
                           rbh_directory_to_distance_matrix,
                           rbh_to_distance_gbgz,
                           ALGrbh_to_algcomboix,
                           plot_umap_from_files,
                           plot_umap_pdf,
                           taxids_to_analyses,
                           taxids_of_interest_to_analyses,
                           topoumap_genmatrix,
                           topoumap_plotumap,
                           sampleToRbhFileDict_to_sample_matrix)

# get the path of this script, so we know where to look for the plotdfs file
# This block imports fasta-parser as fasta
import os
import sys
snakefile_path = os.path.dirname(os.path.realpath(workflow.snakefile))

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

# check that there are some NCBI taxids to plot in the config file
if "taxids" in config:
    config["taxids"] = taxids_to_analyses(config["taxids"])
else:
    config["taxids"] = taxids_of_interest_to_analyses()

config["rbh_directory"] = os.path.abspath(config["rbh_directory"])

# It is possible that this could be called many times. Like once per rule run.
# For this reason this needs to be removed s.t. a list of files is created once.
config["rbh_files"] = list(sorted([os.path.join(config["rbh_directory"], f)
                           for f in os.listdir(config["rbh_directory"])
                           if f.endswith('.rbh')], reverse = True))
print("ACCESSING THE RBH FILES")
config["sample_to_rbh_file"] = {rbh_to_samplename(x, config["ALGname"]): x
                                for x in config["rbh_files"]}

# Results_base_directory is the prefix to which everything will be saved
results_base_directory = "GTUMAP"
if results_base_directory.endswith("/"):
    results_base_directory = results_base_directory[:-1]

# use these parameters for the full space exploration
odog_n    = [10, 15, 20, 35, 50, 75, 100, 150, 250]
odog_m    = [0.0, 0.1, 0.2, 0.5, 0.75, 0.9, 1.0]
odog_size = ["large"]

odol_n    = [5, 15, 50]
odol_m    = [0.0, 0.01, 0.1, 0.2, 0.5, 0.75, 0.9, 1.0]
odol_size = ["large"]
weighting_methods = ["phylogenetic"]

# use these parameters for the clade-specific exploration
codog_n    = [10, 20, 40]
codog_m    = [0.0, 0.25, 0.5, 0.75]
codog_size = ["large"]

rule all:
    input:
        ##    ┓  ┓
        ## ┏┓┏┫┏┓┃ - One-Dot-One-Locus plots
        ## ┗┛┗┻┗┛┗   Each dot represents a single locus, and the data vector is the distance to all other loci
        ## These two sets of files are generated during the rule odolGenCoo and the function topoumap_genmatrix()
        #expand(results_base_directory + "/subchrom/{weighting}/missing_{sizeNaN}/{taxanalysis}.method_{weighting}.missing_{sizeNaN}.coo.npz",
        #        taxanalysis = config["taxids"],
        #        weighting = weighting_methods,
        #        sizeNaN = odol_size),
        #expand(results_base_directory + "/subchrom/{weighting}/missing_{sizeNaN}/{taxanalysis}.method_{weighting}.missing_{sizeNaN}.sampledf.tsv",
        #        taxanalysis = config["taxids"],
        #        weighting = weighting_methods,
        #        sizeNaN = odol_size),
        #expand(results_base_directory + "/subchrom/{weighting}/missing_{sizeNaN}/{taxanalysis}.method_{weighting}.neighbors_{n}.mind_{m}.missing_{sizeNaN}.subchrom.pdf",
        #        n = odol_n,
        #        m = odol_m,
        #        taxanalysis = config["taxids"],
        #        sizeNaN = odol_size,
        #        weighting = weighting_methods),
        #expand(results_base_directory + "/subchrom/{weighting}/missing_{sizeNaN}/{taxanalysis}.method_{weighting}.neighbors_{n}.mind_{m}.missing_{sizeNaN}.subchrom.df",
        #        n = odol_n,
        #        m = odol_m,
        #        taxanalysis = config["taxids"],
        #        sizeNaN = odol_size,
        #        weighting = weighting_methods),
        #expand(results_base_directory + "/subchrom/{weighting}/missing_{sizeNaN}/{taxanalysis}.method_{weighting}.missing_{sizeNaN}.paramsweep.pdf",
        #        taxanalysis = config["taxids"],
        #        sizeNaN = odol_size,
        #        weighting = weighting_methods),
        ##    ┓
        ## ┏┓┏┫┏┓┏┓ - One-Dot-One-Genome plots
        ## ┗┛┗┻┗┛┗┫   Each dot represents a single genome, and the data vector is the distance pairs
        ##        ┛
        #expand(results_base_directory + "/allsamples/allsamples.neighbors_{n}.mind_{m}.missing_{sizeNaN}.df",
        #        n = odog_n,
        #        m = odog_m,
        #        sizeNaN = odog_size),
        #expand(results_base_directory + "/allsamples/allsamples.neighbors_{n}.mind_{m}.missing_{sizeNaN}.bokeh.html",
        #        sample = config["sample_to_rbh_file"].keys(),
        #        n = odog_n,
        #        m = odog_m,
        #        sizeNaN = odog_size),
        #expand(results_base_directory + "/allsamples/allsamples.neighbors_{n}.mind_{m}.missing_{sizeNaN}.pdf",
        #        n = odog_n,
        #        m = odog_m,
        #        sizeNaN = odog_size),
        #expand(results_base_directory + "/allsamples/allsamples.missing_{sizeNaN}.paramsweep.pdf",
        #        sizeNaN = odog_size)
        #    ┓     ┏┓┓ ┏┓┳┓┏┓┏┓
        # ┏┓┏┫┏┓┏┓ ┃ ┃ ┣┫┃┃┣ ┗┓ - One-Dot-One-Genome plots FOR SPECIFIC CLADES
        # ┗┛┗┻┗┛┗┫ ┗┛┗┛┛┗┻┛┗┛┗┛   Each dot represents a single genome, and the data vector is the distance pairs
        #        ┛
        expand(results_base_directory + "/ODOG/clades/{taxanalysis}.neighbors_{n}.mind_{m}.missing_{sizeNaN}.df",
               n       = codog_n,
               m       = codog_m,
               sizeNaN = codog_size,
               taxanalysis = config["taxids"]),
        expand(results_base_directory + "/ODOG/clades/{taxanalysis}.neighbors_{n}.mind_{m}.missing_{sizeNaN}.bokeh.html",
               n       = codog_n,
               m       = codog_m,
               sizeNaN = codog_size,
               taxanalysis = config["taxids"])

# ┏┓    ┓        ┓         ┓
# ┗┓┏┓┏┓┃┏┏┓┏┳┓┏┓┃┏┏┓  ┏┓┓┏┃┏┓┏
# ┗┛┛┗┗┻┛┗┗ ┛┗┗┗┻┛┗┗   ┛ ┗┻┗┗ ┛
#
def generic_get_mem_mb(wildcards, attempt):
    """
    The amount of RAM needed for the script depends on the size of the input genome.
    """
    attemptdict = {1: 2000,
                   2: 4000,
                   3: 8000,
                   4: 16000}
    return attemptdict[attempt]

def pdf_get_mem_mb(wildcards, attempt):
    """
    The amount of RAM needed for the script depends on the size of the input genome.
    """
    attemptdict = {1: 250,
                   2: 1000,
                   3: 4000,
                   4: 10000,
                  }
    return attemptdict[attempt]

def is_valid_gzip(file_path):
    """
    Just check if the gzip file is valid or complete.
    """
    try:
        with gzip.open(file_path, 'rb') as f:
            while f.read(1024 * 1024):  # Read in chunks of 1MB
                pass
        return True
    except (OSError, gzip.BadGzipFile):
        return False

#rule samples_and_gzipped:
#    """
#    This version works in a for loop. I don't like this version because if it fails due to OOM or another reason,
#      it will crash and snakemake will delete all the intermediate files. Instead, it is better that this runs once
#      per file that we need to make.
#
#    The samples need to be converted to a different format before stiching them together into
#      distance matrices in numpy. In an earlier version of this pipeline, I had this rule run once
#      per file. This caused a high IO burden, so now I am running this in a for loop. This actually
#      runs much faster, and each operation takes less than one second typically.
#    """
#    input:
#        rbh_file = [config["sample_to_rbh_file"][sample] for sample in config["sample_to_rbh_file"]]
#    output:
#        gbgz = expand(results_base_directory + "/distance_matrices/{sample}.gb.gz",
#                      sample = config["sample_to_rbh_file"].keys())
#    threads: 1
#    params:
#        ALGname   = config["ALGname"],
#        outprefix = results_base_directory + "/distance_matrices/",
#        basedir   = results_base_directory
#    retries: 4
#    resources:
#        mem_mb = generic_get_mem_mb,
#        highio = 1,
#        runtime = int(len(config["sample_to_rbh_file"]) * 1.5 / 60)
#    run:
#        if not os.path.exists(params.basedir):
#            os.makedirs(params.basedir)
#        if not os.path.exists(params.outprefix):
#            os.makedirs(params.outprefix)
#        num_samples = len(config["sample_to_rbh_file"])
#        counter = 0
#        for thissample in config["sample_to_rbh_file"]:
#            print(f"Reading or converting {thissample} to a distance matrix in gzip format.  {counter}/{num_samples}    ", end = "\r", flush = True)
#            counter += 1
#            inputfile = config["sample_to_rbh_file"][thissample]
#            outfile   = f"{params.outprefix}{thissample}.gb.gz"
#            if os.path.exists(outfile):
#                if is_valid_gzip(outfile):
#                    # If the file is already there and is valid, then we don't need to do anything, and we skip to the next part of the for loop.
#                    continue
#            # print a one-line message to the screen that flushes left each time
#            rbh_to_distance_gbgz(inputfile, outfile, params.ALGname)
#        print()

rule samples_and_gzipped:
    """
    This version works as originally planned. It runs once per file that we need to make.
    If you run too many at once it will freeze the IO on the disk system.

    The samples need to be converted to a different format before stiching them together into
      distance matrices in numpy. In an earlier version of this pipeline, I had this rule run once
      per file. This caused a high IO burden, so now I am running this in a for loop. This actually
      runs much faster, and each operation takes less than one second typically.
    """
    input:
        rbh_file = lambda wildcards: config["sample_to_rbh_file"][wildcards.sample]
    output:
        gbgz = results_base_directory + "/distance_matrices/{sample}.gb.gz"
    threads: 1
    params:
        ALGname   = config["ALGname"]
    retries: 3
    resources:
        mem_mb  = generic_get_mem_mb,
        highio  = 1,
        runtime = 10
    run:
        rbh_to_distance_gbgz(input.rbh_file, output.gbgz, params.ALGname)

rule sample_df_all:
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
        highio = 1,
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
    retries: 3
    resources:
        mem_mb = generic_get_mem_mb,
        highio = 1,
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
                   2: 25000,
                   3: 50000,
                   4: 100000,
                   5: 200000,
                   6: 300000,
                   7: 400000}
    return attemptdict[attempt]

def coo_get_runtime(wildcards, attempt):
    """
    The amount of RAM needed for the script depends on the size of the input genome.
    """
    attemptdict = {1: 60,
                   2: 120,
                   3: 150,
                   4: 180,
                   5: 210,
                   6: 240,
                   7: 270}
    return attemptdict[attempt]

#    ┓
# ┏┓┏┫┏┓┏┓ - One-Dot-One-Genome plots
# ┗┛┗┻┗┛┗┫   Each dot represents a single genome, and the data vector is the distance pairs
#        ┛
rule odogCooGen:
    """
    This generates a coo matrix of the distance matrices for all genomes.
      Each row is a genome, and each column is the distance between every locus pair.
      There is no pre-processing performed on the distance matrices. These are the raw
      distances between every locus pair.

    This takes up a lot of RAM and time. For 3600 genomes, it took:
        CPU Efficiency: 97.32% of 01:57:32 core-walltime
        Job Wall-clock time: 01:57:32
        Memory Utilized: 216.61 GB
        Memory Efficiency: 110.91% of 195.31 GB
    """
    input:
        gbgzfiles = expand(results_base_directory + "/distance_matrices/{sample}.gb.gz",
                           sample = config["sample_to_rbh_file"].keys()),
        combotoindex = results_base_directory + "/combo_to_index.txt"   ,
        sampletsv    = results_base_directory + "/sampledf.tsv"
    output:
        coo    = results_base_directory + "/allsamples.coo.npz"
    threads: 1
    retries: 6
    resources:
        mem_mb = coo_get_mem_mb,
        highio = 10,
        runtime = 300
    run:
        # read in the combo_to_index file as a df, then convert to a dict
        alg_combo_to_ix = algcomboix_file_to_dict(input.combotoindex)
        sampledf = pd.read_csv(input.sampletsv, sep = "\t", index_col = 0)
        coo = construct_coo_matrix_from_sampledf(sampledf, alg_combo_to_ix)
        save_npz(output.coo, coo)

def odogPlotUMAP_get_mem_mb(wildcards, attempt):
    """
    The amount of RAM needed for the script depends on the size of the input genome.
    """
    attemptdict = {#1: 200000,
                   1: 250000,
                   2: 300000,
                   3: 350000,
                   4: 400000}
    return attemptdict[attempt]

def odogPlotUMAP_get_runtime(wildcards, attempt):
    """
    The amount of RAM needed for the script depends on the size of the input genome.
    """
    attemptdict = {#1: 300,
                   1: 360,
                   2: 420,
                   3: 480,
                   4: 540}
    return attemptdict[attempt]

rule odogPlotUMAP:
    input:
        sampletsv    = results_base_directory + "/sampledf.tsv"      ,
        combotoindex = results_base_directory + "/combo_to_index.txt",
        coo    = results_base_directory + "/allsamples.coo.npz"
    output:
        df   = results_base_directory + "/allsamples/allsamples.neighbors_{n}.mind_{m}.missing_{sizeNaN}.df",
        html = results_base_directory + "/allsamples/allsamples.neighbors_{n}.mind_{m}.missing_{sizeNaN}.bokeh.html"
    threads: 1
    params:
        outdir = results_base_directory + "/allsamples",
    retries: 3
    resources:
        mem_mb  = odogPlotUMAP_get_mem_mb,
        runtime = odogPlotUMAP_get_runtime,
        bigUMAPSlots = 1
    run:
        #print(f"These are the wildcards: {wildcards}")
        plot_umap_from_files(input.sampletsv, input.combotoindex, input.coo,
                             params.outdir, "allsamples",
                             wildcards.sizeNaN, int(wildcards.n), float(wildcards.m))

rule odogPDF:
    input:
        df = results_base_directory + "/allsamples/allsamples.neighbors_{n}.mind_{m}.missing_{sizeNaN}.df"
    output:
        pdf = results_base_directory + "/allsamples/allsamples.neighbors_{n}.mind_{m}.missing_{sizeNaN}.pdf"
    threads: 1
    retries: 4
    resources:
        mem_mb = pdf_get_mem_mb,
        highio = 1,
        runtime = 10
    run:
        plot_umap_pdf(input.df, output.pdf, "Allsamples", wildcards.sizeNaN, wildcards.n, wildcards.m)

rule odogSweep:
    """
    Makes a pdf of the parameter sweep of the all-species UMAP plots.
    """
    input:
        dfs = expand(results_base_directory + "/allsamples/allsamples.neighbors_{n}.mind_{m}.missing_{{sizeNaN}}.df",
                     n = odog_n,
                     m = odog_m),
        plotdfs = os.path.join(snakefile_path, "PhyloTreeUMAP_plotdfs.py")
    output:
        pdf = results_base_directory + "/allsamples/allsamples.missing_{sizeNaN}.paramsweep.pdf"
    threads: 1
    retries: 4
    resources:
        mem_mb = pdf_get_mem_mb,
        highio = 1,
        runtime = 5
    shell:
        """
        python {input.plotdfs} -f "{input.dfs}" -o {output.pdf}
        """

#    ┓     ┏┓┓ ┏┓┳┓┏┓┏┓
# ┏┓┏┫┏┓┏┓ ┃ ┃ ┣┫┃┃┣ ┗┓ - One-Dot-One-Genome plots FOR SPECIFIC CLADES
# ┗┛┗┻┗┛┗┫ ┗┛┗┛┛┗┻┛┗┛┗┛   Each dot represents a single genome, and the data vector is the distance pairs
#        ┛
rule odogClCooGen:
    """
    This generates a coo matrix of the distance matrices for the genomes in this clade.

    **This is a modification of the rule called odogCooGen.**
      Each row is a genome, and each column is the distance between every locus pair.
      There is no pre-processing performed on the distance matrices. These are the raw
      distances between every locus pair.

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
        coo       = results_base_directory + "/coo/coo_odog_clade/{taxanalysis}.coo.npz",
        sampletsv = results_base_directory + "/coo/coo_odog_clade/{taxanalysis}.sampledf.tsv"
    threads: 1
    retries: 6
    resources:
        mem_mb = coo_get_mem_mb,
        highio = 1,
        bigUMAPSlots = 1,
        runtime = 300
    params:
        taxids_to_keep    = lambda wildcards: config["taxids"][wildcards.taxanalysis][0],
        taxids_to_remove  = lambda wildcards: config["taxids"][wildcards.taxanalysis][1]
    run:
        # filter the sample df
        sampledf = pd.read_csv(input.sampletsv, sep = "\t", index_col = 0)
        sampledf = filter_sample_df_by_clades(sampledf, params.taxids_to_keep, params.taxids_to_remove)
        # Reset the index. We need to reset the indices so that we can associate the annotated sampledf with the
        #   metadata-less coo matrix.
        sampledf = sampledf.reset_index(inplace = False)
        # save the sampledf. Keep the column.
        sampledf.to_csv(output.sampletsv, sep = "\t")
        print(sampledf)

        # read in the combo_to_index file as a df, then convert to a dict
        alg_combo_to_ix = algcomboix_file_to_dict(input.combotoindex)
        coo = construct_coo_matrix_from_sampledf(sampledf, alg_combo_to_ix)
        save_npz(output.coo, coo)

rule odogClPlotUMAP:
    input:
        sampletsv = results_base_directory + "/coo/coo_odog_clade/{taxanalysis}.sampledf.tsv",
        combotoindex = results_base_directory + "/combo_to_index.txt",
        coo       = results_base_directory + "/coo/coo_odog_clade/{taxanalysis}.coo.npz",
    output:
        df   = results_base_directory + "/ODOG/clades/{taxanalysis}.neighbors_{n}.mind_{m}.missing_{sizeNaN}.df",
        html = results_base_directory + "/ODOG/clades/{taxanalysis}.neighbors_{n}.mind_{m}.missing_{sizeNaN}.bokeh.html"
    threads: 1
    retries: 3
    resources:
        mem_mb  = odogPlotUMAP_get_mem_mb,
        runtime = odogPlotUMAP_get_runtime,
        bigUMAPSlots = 1
    run:
        #print(f"These are the wildcards: {wildcards}")
        plot_umap_from_files(input.sampletsv, input.combotoindex, input.coo,
                             wildcards.taxanalysis, wildcards.sizeNaN, int(wildcards.n),
                             float(wildcards.m), output.df, output.html)

#    ┓  ┓
# ┏┓┏┫┏┓┃ - One-Dot-One-Locus plots
# ┗┛┗┻┗┛┗   Each dot represents a single locus, and the data vector is the distance to all other loci
#
rule odolGenCoo:
    input:
        sampletsv    = results_base_directory + "/sampledf.tsv",
        combotoindex = results_base_directory + "/combo_to_index.txt",
        coo          = results_base_directory + "/allsamples.coo.npz",
        ALGrbh       = config["ALG_rbh_file"]
    output:
        coo      = results_base_directory + "/subchrom/{weighting}/missing_{sizeNaN}/{taxanalysis}.method_{weighting}.missing_{sizeNaN}.coo.npz",
        sampledf = results_base_directory + "/subchrom/{weighting}/missing_{sizeNaN}/{taxanalysis}.method_{weighting}.missing_{sizeNaN}.sampledf.tsv"
    threads: 1
    params:
        outdir = results_base_directory + "/allsamples",
        taxids_to_keep    = lambda wildcards: config["taxids"][wildcards.taxanalysis][0],
        taxids_to_remove  = lambda wildcards: config["taxids"][wildcards.taxanalysis][1]
    retries: 5
    resources:
        mem_mb = coo_get_mem_mb,
        highio = 1,
        runtime = 60
    run:
        topoumap_genmatrix(input.sampletsv, input.combotoindex, input.coo, input.ALGrbh,
                           wildcards.taxanalysis, params.taxids_to_keep, params.taxids_to_remove,
                           output.coo, output.sampledf,
                           wildcards.sizeNaN,
                           method = wildcards.weighting,
                           missing_value_as = 9999999999)

def odol_plot_get_mem_mb(wildcards, attempt):
    """
    The amount of RAM needed for the script depends on the size of the input genome.
    """
    attemptdict = {1: 2000,
                   2: 30000,
                   3: 50000,
                   4: 100000,
                   5: 200000,
                   6: 300000,
                   7: 400000}
    return attemptdict[attempt]

def odol_plot_get_runtime(wildcards, attempt):
    """
    The amount of RAM needed for the script depends on the size of the input genome.
    """
    attemptdict = {1: 5,
                   2: 120,
                   3: 150,
                   4: 180,
                   5: 210,
                   6: 240,
                   7: 270}
    return attemptdict[attempt]

rule odolPlotUMAP:
    input:
        sampletsv    = results_base_directory + "/sampledf.tsv",
        coo = results_base_directory + "/subchrom/{weighting}/missing_{sizeNaN}/{taxanalysis}.method_{weighting}.missing_{sizeNaN}.coo.npz",
        ALGrbh = config["ALG_rbh_file"]
    output:
        df     = results_base_directory + "/subchrom/{weighting}/missing_{sizeNaN}/{taxanalysis}.method_{weighting}.neighbors_{n}.mind_{m}.missing_{sizeNaN}.subchrom.df",
        html   = results_base_directory + "/subchrom/{weighting}/missing_{sizeNaN}/{taxanalysis}.method_{weighting}.neighbors_{n}.mind_{m}.missing_{sizeNaN}.subchrom.bokeh.html",
        jpeg   = results_base_directory + "/subchrom/{weighting}/mixxing_{sizeNaN}/{taxanalysis}.method_{weighting}.neighbors_{n}.mind_{m}.missing_{sizeNaN}.subchrom.connectivity.jpeg"
    params:
        outdir = results_base_directory + "/subchrom",
    retries: 6
    threads: 1
    resources:
        mem_mb = odol_plot_get_mem_mb,
        highio = 1,
        runtime = odol_plot_get_runtime
    run:
        topoumap_plotumap(wildcards.taxanalysis, input.sampletsv, input.ALGrbh, input.coo,
                          params.outdir, wildcards.sizeNaN, int(wildcards.n), float(wildcards.m),
                          output.df, output.html, output.jpeg)

def odolPlotPdf_get_mem_mb(wildcards, attempt):
    """
    The amount of RAM needed for the script depends on the size of the input genome.
    """
    attemptdict = {1: 1000,
                   2: 2000,
                   3: 4000,
                   4: 8000,
                   5: 16000,
                  }
    return attemptdict[attempt]

rule odolPlotPdf:
    input:
        df  = results_base_directory + "/subchrom/{weighting}/missing_{sizeNaN}/{taxanalysis}.method_{weighting}.neighbors_{n}.mind_{m}.missing_{sizeNaN}.subchrom.df",
    output:
        pdf = results_base_directory + "/subchrom/{weighting}/missing_{sizeNaN}/{taxanalysis}.method_{weighting}.neighbors_{n}.mind_{m}.missing_{sizeNaN}.subchrom.pdf"
    threads: 1
    retries: 4
    resources:
        mem_mb = odolPlotPdf_get_mem_mb,
        highio = 1,
        runtime = 20
    run:
        plot_umap_pdf(input.df, output.pdf, wildcards.taxanalysis, wildcards.sizeNaN, wildcards.n, wildcards.m)

rule odolSweep:
    """
    Makes a pdf of the parameter sweep of the clade-specific UMAP plots.
    """
    input:
        dfs = expand(results_base_directory + "/subchrom/{{weighting}}/missing_{{sizeNaN}}/{{taxanalysis}}.method_{{weighting}}.neighbors_{n}.mind_{m}.missing_{{sizeNaN}}.subchrom.df",
                      n = odol_n,
                      m = odol_m),
        plotdfs = os.path.join(snakefile_path, "PhyloTreeUMAP_plotdfs.py")
    output:
        pdf = results_base_directory + "/subchrom/{weighting}/missing_{sizeNaN}/{taxanalysis}.method_{weighting}.missing_{sizeNaN}.paramsweep.pdf"
    threads: 1
    retries: 4
    resources:
        mem_mb = pdf_get_mem_mb,
        highio = 1,
        runtime = 5
    shell:
        """
        python {input.plotdfs} -f "{input.dfs}" -o {output.pdf}
        """