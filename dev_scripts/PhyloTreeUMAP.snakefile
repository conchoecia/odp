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
                           plot_umap_pdf,
                           topoumap_genmatrix,
                           topoumap_plotumap,
                           sampleToRbhFileDict_to_sample_matrix)

# get the path of this script, so we know where to look for the plotdfs file
# This block imports fasta-parser as fasta
import os
import sys
snakefile_path = os.path.dirname(os.path.realpath(workflow.snakefile))

from ete3 import NCBITaxa,Tree
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

# check that there are some NCBI taxids to plot in the config file
if "taxids" not in config:
    # 33317 is protostomes
    # 33213 is bilateria
    config["taxids"] = [ [[10197], []],      # ctenophores
                         [[6040],  [60882]], # porifera minus Hexactinellida
                         [[6073],  []],      # cnidaria
                         [[6340],  [42113]], # annelida minus clitellata
                         [[42113], []],     # clitellata
                         [[6447],  [6606]],  # mollusca minus coleoida
                         [[6606],  []],      # coleoida
                         [[50557], []],     # insecta
                         [[32341], []],     # Sophophora - subset of drosophilids
                         #[[61985], []],     # myriapoda
                         [[6231],  []],      # nematoda
                         [[7586],  []],      # echinodermata
                         [[7742],  []],      # Vertebrata
                         #[[33317],[]]
                        ]
# Come up with the taxid analyses. Each entry will have a string indicating what is in it and what is not.
# Bilateria_33213_without_None if we want to plot all bilateria, and want to remove specific things
# Bilateria_33213_without_33317_7652 if we want to plot all bilateria, but we don't want to plot the protostomes or lytechinus
# Bilateria_33213_without_33317_7652 if we want to plot all bilateria, but we don't want to plot the protostomes or lytechinus
analyses = {}
ncbi = NCBITaxa()
for entry in config["taxids"]:
    # get the clade name to make reading easier
    clade = ncbi.get_taxid_translator([entry[0][0]])[entry[0][0]].replace(" ", "").replace("-", "").replace(".", "")
    # make sure that the length of the 0th entry is at least length 1
    if len(entry[0]) == 0:
        raise IOError("There must be at least one taxid in the first entry of the taxids list.")
    analysis_name = clade + "_" + "_".join([str(x) for x in entry[0]]) + "_without_"
    analysis_name += "_".join([str(x) for x in entry[1]]) if len(entry[1]) > 0 else "None"
    analyses[analysis_name] = entry
config["taxids"] = analyses
print(config["taxids"])

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

odog_n    = [10, 15, 20, 25, 30, 40, 50, 75, 100, 150, 250]
odog_m    = [0.0, 0.1, 0.2, 0.5, 0.65, 0.75, 0.9, 1.0]
odog_size = ["small"]

odol_n    = [5, 10, 15, 50]
odol_m    = [0.0, 0.01, 0.1, 0.2, 0.5, 0.75, 0.9, 1.0]
odol_size = ["small"]
rule all:
    input:
        #    ┓  ┓
        # ┏┓┏┫┏┓┃ - One-Dot-One-Locus plots
        # ┗┛┗┻┗┛┗   Each dot represents a single locus, and the data vector is the distance to all other loci
        #
        expand(results_base_directory + "/subchrom/{taxanalysis}.coo.npz",
                taxanalysis = config["taxids"]),
        expand(results_base_directory + "/subchrom/{taxanalysis}.neighbors_{n}.mind_{m}.missing_{sizeNaN}.subchrom.pdf",
                n = odol_n,
                m = odol_m,
                taxanalysis = config["taxids"],
                sizeNaN = odol_size),
        expand(results_base_directory + "/subchrom/{taxanalysis}.neighbors_{n}.mind_{m}.missing_{sizeNaN}.subchrom.df",
                n = odol_n,
                m = odol_m,
                taxanalysis = config["taxids"],
                sizeNaN = odol_size),
        expand(results_base_directory + "/subchrom/{taxanalysis}.missing_{sizeNaN}.paramsweep.pdf",
                taxanalysis = config["taxids"],
                sizeNaN = odol_size),
        #    ┓
        # ┏┓┏┫┏┓┏┓ - One-Dot-One-Genome plots
        # ┗┛┗┻┗┛┗┫   Each dot represents a single genome, and the data vector is the distance pairs
        #        ┛
        expand(results_base_directory + "/allsamples/allsamples.neighbors_{n}.mind_{m}.missing_{sizeNaN}.df",
                n = odog_n,
                m = odog_m,
                sizeNaN = odog_size),
        expand(results_base_directory + "/allsamples/allsamples.neighbors_{n}.mind_{m}.missing_{sizeNaN}.bokeh.html",
                sample = config["sample_to_rbh_file"].keys(),
                n = odog_n,
                m = odog_m,
                sizeNaN = odog_size),
        expand(results_base_directory + "/allsamples/allsamples.neighbors_{n}.mind_{m}.missing_{sizeNaN}.pdf",
                n = odog_n,
                m = odog_m,
                sizeNaN = odog_size),
        expand(results_base_directory + "/allsamples/allsamples.missing_{sizeNaN}.paramsweep.pdf",
                sizeNaN = odog_size)


# ┏┓    ┓        ┓         ┓
# ┗┓┏┓┏┓┃┏┏┓┏┳┓┏┓┃┏┏┓  ┏┓┓┏┃┏┓┏
# ┗┛┛┗┗┻┛┗┗ ┛┗┗┗┻┛┗┗   ┛ ┗┻┗┗ ┛
#
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
                   2: 25000,
                   3: 50000,
                   4: 100000,
                   5: 200000,
                   6: 300000,
                  }
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
                  }
    return attemptdict[attempt]


#    ┓
# ┏┓┏┫┏┓┏┓ - One-Dot-One-Genome plots
# ┗┛┗┻┗┛┗┫   Each dot represents a single genome, and the data vector is the distance pairs
#        ┛
rule odogCooGen:
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
    retries: 6
    resources:
        mem_mb = coo_get_mem_mb,
        runtime = 300
    run:
        # read in the combo_to_index file as a df, then convert to a dict
        alg_combo_to_ix = algcomboix_file_to_dict(input.combotoindex)
        sampledf = pd.read_csv(input.sampletsv, sep = "\t", index_col = 0)
        coo = construct_coo_matrix_from_sampledf(sampledf, alg_combo_to_ix)
        save_npz(output.coo, coo)

rule odogPlotUMAP:
    input:
        sampletsv    = results_base_directory + "/sampledf.tsv",
        combotoindex = results_base_directory + "/combo_to_index.txt",
        coo    = results_base_directory + "/allsamples.coo.npz"
    output:
        df   = results_base_directory + "/allsamples/allsamples.neighbors_{n}.mind_{m}.missing_{sizeNaN}.df",
        html = results_base_directory + "/allsamples/allsamples.neighbors_{n}.mind_{m}.missing_{sizeNaN}.bokeh.html"
    threads: 1
    params:
        outdir = results_base_directory + "/allsamples",
    retries: 6
    resources:
        mem_mb = coo_get_mem_mb,
        runtime = coo_get_runtime
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
        runtime = 5
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
        runtime = 5
    shell:
        """
        python {input.plotdfs} -f "{input.dfs}" -o {output.pdf}
        """

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
        coo = results_base_directory + "/subchrom/{taxanalysis}.coo.npz"
    threads: 1
    params:
        outdir = results_base_directory + "/allsamples",
        taxids_to_keep    = lambda wildcards: config["taxids"][wildcards.taxanalysis][0],
        taxids_to_remove  = lambda wildcards: config["taxids"][wildcards.taxanalysis][1]
    retries: 5
    resources:
        mem_mb = coo_get_mem_mb,
        runtime = 60
    run:
        topoumap_genmatrix(input.sampletsv, input.combotoindex, input.coo, input.ALGrbh,
                           wildcards.taxanalysis, params.taxids_to_keep, params.taxids_to_remove,
                           output.coo)

rule odolPlotUMAP:
    input:
        sampletsv    = results_base_directory + "/sampledf.tsv",
        coo = results_base_directory + "/subchrom/{taxanalysis}.coo.npz",
        ALGrbh = config["ALG_rbh_file"]
    output:
        df     = results_base_directory + "/subchrom/{taxanalysis}.neighbors_{n}.mind_{m}.missing_{sizeNaN}.subchrom.df",
        html   = results_base_directory + "/subchrom/{taxanalysis}.neighbors_{n}.mind_{m}.missing_{sizeNaN}.subchrom.bokeh.html"
    params:
        outdir = results_base_directory + "/subchrom",
    retries: 5
    threads: 1
    resources:
        mem_mb = coo_get_mem_mb,
        runtime = 120
    run:
        topoumap_plotumap(wildcards.taxanalysis, input.sampletsv, input.ALGrbh, input.coo,
                          params.outdir, wildcards.sizeNaN, int(wildcards.n), float(wildcards.m))

rule odolPlotPdf:
    input:
        df  = results_base_directory + "/subchrom/{taxanalysis}.neighbors_{n}.mind_{m}.missing_{sizeNaN}.subchrom.df",
    output:
        pdf = results_base_directory + "/subchrom/{taxanalysis}.neighbors_{n}.mind_{m}.missing_{sizeNaN}.subchrom.pdf"
    threads: 1
    retries: 4
    resources:
        mem_mb = pdf_get_mem_mb,
        runtime = 5
    run:
        plot_umap_pdf(input.df, output.pdf, wildcards.taxanalysis, wildcards.sizeNaN, wildcards.n, wildcards.m)

rule odolSweep:
    """
    Makes a pdf of the parameter sweep of the clade-specific UMAP plots.
    """
    input:
        dfs = expand(results_base_directory + "/subchrom/{{taxanalysis}}.neighbors_{n}.mind_{m}.missing_{{sizeNaN}}.subchrom.df",
                      n = odol_n,
                      m = odol_m),
        plotdfs = os.path.join(snakefile_path, "PhyloTreeUMAP_plotdfs.py")
    output:
        pdf = results_base_directory + "/subchrom/{taxanalysis}.missing_{sizeNaN}.paramsweep.pdf"
    threads: 1
    retries: 4
    resources:
        mem_mb = pdf_get_mem_mb,
        runtime = 5
    shell:
        """
        python {input.plotdfs} -f "{input.dfs}" -o {output.pdf}
        """