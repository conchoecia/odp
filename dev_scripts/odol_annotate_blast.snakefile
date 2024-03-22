"""
Program  : odol_annotate_blast.snakefile
Language : snakemake
Date     : 2024-03-06
Author   : Darrin T. Schultz
Email    : darrin.schultz@univie.ac.at
Github   : https://github.com/conchoecia/odp
Support  : For issues or questions, please search if the topic has been discussed already
           on github and open a new issue if not: https://github.com/conchoecia/odp/issues
License  : GNU GENERAL PUBLIC LICENSE, Version 3, 29 June 2007. See the LICENSE file.
Citation : No citation currently.

Description:
  - Takes in the subchrom directory of PhyloTreeUMAP.snakefile
    - The directory contains the *.coo.npz files that are the sparse matrices of each set of species.
    - The directory also has the *.sampledf.tsv files that have the info about all of the species in the set.
      - This df will be used to get the samples for each set of species.
      - The samples for each set of species must be present in the config file, as we will use their proteins.

Usage instructions:
  - None currently
"""

# import the combinations
from itertools import product
# plotting stuff
import bokeh
# needs to be umap v0.5.5 or later
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import umap
import umap.plot
import random
import time
import warnings

from ast import literal_eval
import os
import pandas as pd
from scipy.sparse import coo_matrix, lil_matrix, save_npz, load_npz, csr_matrix
import sys

from PhyloTreeUMAP import (algcomboix_file_to_dict,
                           construct_coo_matrix_from_sampledf,
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

from odol_annotate_blast import umapdf_one_species_one_query

# ODP-specific imports
snakefile_path = os.path.dirname(os.path.realpath(workflow.snakefile))
scripts_path = os.path.join(snakefile_path, "../scripts")
dependencies_path = os.path.join(snakefile_path, "../dependencies")
sys.path.insert(1, scripts_path)
import odp_functions as odpf
from rbh_tools import parse_rbh

configfile: "config.yaml"
basedir = "UMAP_blast_results"
# First we define all of the RBH files
# get the rbh files in the directory
config["ALGname"] = "BCnSSimakov2022"

# CHECKS
# CHECK SUBCHROM_DIRECTORY - Check that rbh_directory exists in the config file
if not "subchrom_directory" in config:
    raise ValueError("subchrom_directory not in config file")
# Make sure that the directory exists
if not os.path.exists(config["subchrom_directory"]):
    raise ValueError(f"subchrom_directory {config['subchrom_directory']} does not exist")

# CHECK RBH FILE - Make sure that the user has provided the rbh file
if not "ALG_rbh_file" in config:
    raise ValueError("ALG_rbh_file not in config file")
# Check that the ALG rbh file exists
if not os.path.exists(config["ALG_rbh_file"]):
    raise ValueError(f"rbh_directory {config['rbh_directory']} does not exist")

# CHECK DUPLICATE PROTEINS - Make sure that the "duplicate_proteins" is in the config file
if not "duplicate_proteins" in config:
    raise ValueError("duplicate_proteins not in config file")

# CHECK blast_files
if not "blast_files" in config:
    raise ValueError("blast_files not in config file")
for x in config["blast_files"]:
    if not os.path.exists(config['blast_files'][x]):
        raise ValueError(f"blast_files {config['blast_files'][x]} does not exist")

# Check first colon in fasta name
if "blast_split_first_colon" not in config:
    config["blast_split_first_colon"] = False
else:
    if not isinstance(config["blast_split_first_colon"], bool):
        raise ValueError("blast_split_first_colon must be a boolean")

# We need to define config["species_of_interest"] in case we want to make any plots for single species.
# We probably don't want to make thousands of pdfs of all the clusters to look at, so just do this for some species we care about at the moment.
# For now, do this by defining the species of interest in the config file.
# If the "species_of_interest" is not in the config file, then we initialize it as empty.
if "species_of_interest" not in config:
    config["species_of_interest"] = []

# Make an analysis_to_sampledf and analysis_to_coo dictionary.
#  - The keys are the analysis names.
#  - The values are the filepath to the sampledf file.
coo_files            = [x for x in os.listdir(config["subchrom_directory"]) if x.endswith(".coo.npz")]
if len(coo_files) == 0:
    raise ValueError(f"No .coo.npz files in {config['subchrom_directory']}")
analysis_to_coo      = {x.split(".method_")[0].replace(".coo.npz", ""): os.path.join(config["subchrom_directory"], x)
                        for x in coo_files}
sample_files = [x for x in os.listdir(config["subchrom_directory"]) if x.endswith(".sampledf.tsv")]
analysis_to_sampledf = {}
for thisanalysis in analysis_to_coo:
    thissampledf = [x for x in sample_files if thisanalysis in x][0]
    analysis_to_sampledf[thisanalysis] = os.path.join(config["subchrom_directory"], thissampledf)

# Load in all of the dfs and make an analysis_to_samples dictionary.
#  - The keys are the analysis names.
#  - The values are lists of sample names that belong in that analysis
analysis_to_samples = {k: list(pd.read_csv(v, sep="\t", index_col=0)["sample"])
                       for k, v in analysis_to_sampledf.items()}
unique_samples = list(set([x for y in analysis_to_samples.values() for x in y]))
sample_to_analysis = {}
for thisanalysis in analysis_to_samples:
    for thissample in analysis_to_samples[thisanalysis]:
        sample_to_analysis[thissample] = thisanalysis

for k in analysis_to_sampledf:
    print(f"{k}: {len(analysis_to_samples[k])} genomes")

# Check RBH files for each sample
if "rbh_directory" not in config:
    raise ValueError("rbh_directory not in config file")
if not os.path.exists(config["rbh_directory"]):
    raise ValueError(f"rbh_directory {config['rbh_directory']} does not exist")
# traverse the files in the directory that end in .rbh and make a sample_to_rbh_file dictionary
rbh_files = [x for x in os.listdir(config["rbh_directory"]) if x.endswith(".rbh")]
sample_to_rbh_file = {x.split("_")[1]: os.path.join(config["rbh_directory"], x) for x in rbh_files}
# check that all of the unique_samples are in the sample_to_rbh_file directory
if not all([x in sample_to_rbh_file for x in unique_samples]):
    raise ValueError("Not all of the unique_samples are in the rbh_directory")

#print(analysis_to_samples)
#print(sample_to_analysis)

#odol_n = [5, 10, 15]
#odol_m = [0, 0.01, 0.1]
odol_n = [15]
odol_m = [0.1]
# we need to make a list of target files, because expand doesn't have the ability for us
# to subset the analyses we need for plotting the one-analysis-one-query-one-sample plots
file_targets = []
for thissp in config["species_of_interest"]:
    if thissp in sample_to_analysis:
        thisanalysis = sample_to_analysis[thissp]
        # Now we have handled the species, and the analysis that needs to be plotted.
        # We need to get all the unique combinations of the odol_n, odol_m, ALG, and query
        for n, m, ALG, query in product(odol_n, odol_m, config["targetALGs"], config["blast_files"]):
            f = basedir + f"/ALG_reimbedding/one_analysis_one_species_one_query/{thisanalysis}.{thissp}.{query}.{ALG}.neighbors_{n}.mind_{m}.missing_large.subchrom.pdf"
            file_targets.append(f)
print("species of interest are: {}".format(config["species_of_interest"]))
print("analysis_to_sampledf: {}".format(analysis_to_sampledf))

wildcard_constraints:
    analysis = "[A-Za-z0-9_]+",
    query    = "[A-Za-z0-9]+" ,
    ALG      = "[A-Za-z0-9_]+",

rule all:
    input:
        expand(basedir + "/blast_concat/{query}/{analysis}.concat.filt.blastp",
                query = config["blast_files"],
                analysis = analysis_to_samples.keys()),
        expand(basedir + "/filt_coo/{analysis}_{ALG}.matrix.filt.tsv.gz",
                analysis = analysis_to_samples.keys(),
                ALG = config["targetALGs"]),
        #expand(basedir + "/ALG_reimbedding/{analysis}/{ALG}/{analysis}.neighbors_{n}.mind_{m}.missing_large.subchrom.df",
        #        analysis = analysis_to_samples.keys(),
        #        ALG = config["targetALGs"],
        #        n = odol_n,
        #        m = odol_m),
        #expand(basedir + "/ALG_reimbedding/{analysis}/{ALG}/{analysis}.neighbors_{n}.mind_{m}.missing_large.subchrom.bokeh.html",
        #        analysis = analysis_to_samples.keys(),
        #        ALG = config["targetALGs"],
        #        n = odol_n,
        #        m = odol_m),
        ## pre-plotting df
        #expand(basedir + "/ALG_reimbedding/{analysis}/{ALG}/{analysis}.neighbors_{n}.mind_{m}.missing_large.subchrom.{query}.df",
        #        query = config["blast_files"],
        #        analysis = analysis_to_samples.keys(),
        #        ALG = config["targetALGs"],
        #        n = odol_n,
        #        m = odol_m),
        #expand(basedir + "/ALG_reimbedding/{analysis}/{ALG}/{analysis}.neighbors_{n}.mind_{m}.missing_large.subchrom.{query}.pdf",
        #        query = config["blast_files"],
        #        analysis = analysis_to_samples.keys(),
        #        ALG = config["targetALGs"],
        #        n = odol_n,
        #        m = odol_m),
        # file_targets contains the one-analysis-one-query-one-sample plots
        file_targets

rule install_diamond:
    output:
        diamond = os.path.join(dependencies_path, "diamond/diamond"),
        tar = temp("diamond-linux64.tar.gz")
    threads: 1
    resources:
        mem_mb  = 1000,
        runtime = 3
    shell:
        """
        wget http://github.com/bbuchfink/diamond/releases/download/v2.1.9/diamond-linux64.tar.gz
        tar xzf diamond-linux64.tar.gz
        mv diamond {output.diamond}
        """

def legality_get_mem_mb(wildcards, attempt):
    """
    The amount of RAM needed could change.
    """
    attemptdict = {1: 4000,
                   2: 8000,
                   3: 16000,
                   4: 64000,
                   5: 256000,
                   6: 512000}
    return attemptdict[attempt]

rule check_input_legality:
    """
    The very first thing that we must do is to check the input files.
    Determine whether the input files are formatted correctly.
    We call a function that checks these files:
      - The genome .fasta file
      - The .chrom file
      - The protein .fasta file

    Notes:
    202140106 - This job appears to request too much RAM, causing OOM errors even
        when the input files are small. I'm not sure why this is happening.
        As a result, I am now adding a resources to limit the number of simultaneous jobs.
    """
    input:
        fasta = lambda wildcards: config["species"][wildcards.sample]["genome"],
        chrom = lambda wildcards: config["species"][wildcards.sample]["chrom"],
        pep   = lambda wildcards: config["species"][wildcards.sample]["proteins"]
    output:
        input_pass = basedir + "/db/input_check/{sample}_pass.txt"
    params:
        duplicate_handling = config["duplicate_proteins"]
    resources:
        mem_mb = legality_get_mem_mb,
        runtime   = 30
    threads: 1
    retries: 6
    run:
        # We always check to see if the chrom file is legal, regardless of whether we're ignoring duplicate prots.
        if not odpf.chrom_file_is_legal(input.chrom):
            raise IOError("The chrom file is not legal: {}".format(input.chrom))
        if params.duplicate_handling == "fail":
            # We want to enforce that the input files are legal
            # See the odpf.check_species_input_legality function for more details
            if odpf.check_species_input_legality(input.fasta, input.pep, input.chrom):
                with open(output.input_pass, "w") as outf:
                    outf.write("pass")
        elif params.duplicate_handling == "pass":
            # we just assume that everything is copacetic and let the program continue
            with open(output.input_pass, "w") as outf:
                outf.write("pass")
        else:
            raise IOError("Right now we can only handle fail and pass cases. Sorry.")

rule make_diamond_and_blast_db:
    """
    We make both diamond and blastp databases for the proteins.
    We only do this one the sample has passed the input QC checks.
    """
    input:
        pep = lambda wildcards: config["species"][wildcards.sample]["proteins"],
        input_pass = basedir + "/db/input_check/{sample}_pass.txt",
        diamond = os.path.join(dependencies_path, "diamond/diamond")
    output:
        dmnd       = basedir + "/db/dmnd/{sample}_prots.gz.dmnd"
    params:
        outdir = basedir + "/db/dmnd/"
    resources:
        mem_mb = 1000, # should be small
        runtime   = 10  # ten minutes
    threads: 1 # No reason to have more than one thread
    shell:
        """
        mkdir -p {params.outdir}
        {input.diamond} makedb --in {input.pep} --db {output.dmnd}
        """

rule blast_against_sample:
    input:
        query = lambda wildcards: config["blast_files"][wildcards.query],
        dmnd  = basedir + "/db/dmnd/{sample}_prots.gz.dmnd",
        diamond = os.path.join(dependencies_path, "diamond/diamond")
    output:
        blastp  = basedir + "/blast/{query}/{sample}_results.blastp"
    threads: 8
    resources:
        mem_mb = 8000,
        runtime = 10
    shell:
        """
        {input.diamond} blastp --query {input.query} --db {input.dmnd} --threads {threads} --evalue 1E-5 --outfmt 6 --out {output.blastp}
        """

def tophit_get_mem_mb(wildcards, attempt):
    """
    The amount of RAM needed could change.
    """
    attemptdict = {1: 250,
                   2: 500,
                   3: 1000,
                   4: 2000,
                   5: 4000,
                   6: 8000}
    return attemptdict[attempt]

rule top_hit:
    input:
        blastp       = basedir + "/blast/{query}/{sample}_results.blastp",
        rbh_file     = lambda wildcards: sample_to_rbh_file[wildcards.sample],
        alg_rbh_file = config["ALG_rbh_file"],
        sample_chrom = lambda wildcards: config["species"][wildcards.sample]["chrom"]
    output:
        blastp  = basedir + "/blast_filt/{query}/{sample}_results.filt.blastp"
    retries: 6
    threads: 1
    resources:
        mem_mb = tophit_get_mem_mb,
        runtime = 1
    params:
        split_first_colon = config["blast_split_first_colon"]
    run:
        # parse the ALG rbh file
        rbh_df = parse_rbh(input.alg_rbh_file)
        genefamilies = list(rbh_df["rbh"])
        # check if input.blastp is empty
        if os.stat(input.blastp).st_size == 0:
            # make an empty output.blastp file
            with open(output.blastp, "w") as f:
                pass
        else:
            # use pandas and remove rows that are duplicates of the first column
            df = pd.read_csv(input.blastp, sep="\t", header=None)
            df = df.drop_duplicates(subset=0)
            # correct an error where, if there are tab characters in the protein database, the tab characters are replaced in the protein names
            df[0] = df[0].apply(lambda x: x.split("\\t")[0])

            # if we split on the first colon, do that not
            if params.split_first_colon:
                df[0] = df[0].apply(lambda x: x.split(":")[0])
                # now sort by column 0, then column 11 (ascending, descending)
                df = df.sort_values(by=[0, 11], ascending=[True, False])
                # remove duplicates that we see from the first column
                df = df.drop_duplicates(subset=0)
            # Now we check if these were annotated with rbh or not.
            # If they were, we delete all of those entries.
            db_rbhs = set(["_".join(x.split("_")[:-1]) for x in df[1]])
            # If more than 25% of the things in db_rbhs are in genefamilies, this was annotated with rbh
            #  and we delete everything.
            keepgoing = True
            if len(db_rbhs.intersection(genefamilies)) / len(db_rbhs) > 0.25:
                # delete all rows
                df = df.drop(df.index)
                keepgoing = False

            if keepgoing:
                # now we add some extra columns. They are: ["scaf_of_nearest_ALG", "nearest_ALGs", "position_of_nearest_ALG", "dist_to_nearest_ALG"]
                #  - "scaf_of_nearest_ALG"     - this is a single value
                #  - "nearest_ALGs"            - a list of the closest to furthest away orthologs that exist on th same chromosomes
                #  - "position_of_nearest_ALG" - a list of the positions of the closest to furthest away orthologs that exist on the same chromosomes
                #  - "dist_to_nearest_ALG"     - a list of the distances to the closest to furthest away orthologs that exist on the same chromosomes
                sample_rbh = parse_rbh(input.rbh_file)
                sample_chrom = pd.read_csv(input.sample_chrom, sep="\t")
                sample_chrom.columns = ["protein", "scaf", "strand", "start", "stop"]
                df["scaf_of_nearest_ALG"]     = "-1"
                df["nearest_ALG"]             = [[] for x in range(len(df))]
                df["position_of_nearest_ALG"] = [[] for x in range(len(df))]
                df["dist_to_nearest_ALG"]     = [[] for x in range(len(df))]
                for i, row in df.iterrows():
                    thisprot = row[1]
                    nullrow = False
                    if thisprot in sample_chrom["protein"].values:
                        # get the scaffold and position of the db protein
                        #print(f"The sample is {wildcards.sample} the prot is {row[1]}", file = sys.stderr)
                        prot_scaf = sample_chrom.loc[sample_chrom["protein"] == thisprot, "scaf"].values[0]
                        # get the position of the db protein
                        prot_pos = sample_chrom.loc[sample_chrom["protein"] == thisprot, "start"].values[0]

                        # subset the rbh df to just have entries that are on the same scaffold as this protein
                        subrbh = sample_rbh[sample_rbh[f"{wildcards.sample}_scaf"] == prot_scaf].copy()
                        subrbh["dist_to_prot"] = abs(subrbh[f"{wildcards.sample}_pos"] - prot_pos)
                        # sort by the distance to the protein
                        subrbh = subrbh.sort_values(by="dist_to_prot", ascending=True)
                        if len(subrbh) > 0:
                            df.at[i, "scaf_of_nearest_ALG"] = prot_scaf
                            df.at[i, "nearest_ALG"] = list(subrbh["rbh"])
                            df.at[i, "dist_to_nearest_ALG"] = list(subrbh["dist_to_prot"])
                            df.at[i, "position_of_nearest_ALG"] = list(subrbh[f"{wildcards.sample}_pos"])
            df.to_csv(output.blastp, sep="\t", index=False, header=False)

rule collateAnalysis:
    input:
        blast_entries  = lambda wildcards: expand(basedir + "/blast_filt/{{query}}/{sample}_results.filt.blastp",
                                sample = analysis_to_samples[wildcards.analysis])
    output:
        concatenated = basedir + "/blast_concat/{query}/{analysis}.concat.filt.blastp"
    params:
        blastprefix = basedir + "/blast_filt/{query}/"
    threads: 1
    resources:
        mem_mb = tophit_get_mem_mb,
        runtime = 1
    retries: 6
    run:
        print (f"The analysis is {wildcards.analysis}")
        print(f"The samples are {analysis_to_samples[wildcards.analysis]}")
        # concatenate all of the entries and add the sample name to the first column
        with open(output.concatenated, "w") as outf:
            for sample in analysis_to_samples[wildcards.analysis]:
                print(f"This sample is {sample}")
                blastfile = os.path.join(params.blastprefix, f"{sample}_results.filt.blastp")
                # add the sample name to the first column
                with open(blastfile, "r") as f:
                    for line in f:
                        line = line.strip()
                        if line:
                            splitd = [sample] + line.split("\t")
                            splitd = [x.split("\\t")[0] for x in splitd if x != ""]
                            print("\t".join(splitd), file=outf)

rule filtCOO:
    """
    For every clade, takes the coo file and filters it based on which ALGs we want to keep.
    """
    input:
        coo = lambda wildcards: analysis_to_coo[wildcards.analysis],
        alg_rbh_file = config["ALG_rbh_file"]
    output:
        filt_coo = basedir + "/filt_coo/{analysis}_{ALG}.matrix.filt.tsv.gz"
    threads: 1
    params:
        thisALG = lambda wildcards: wildcards.ALG
    resources:
        mem_mb = tophit_get_mem_mb,
        runtime = 1
    retries: 6
    run:
        # readh in the alg rbh dataframe
        algrbhdf = parse_rbh(input.alg_rbh_file)
        # Get the rows where the gene_group column matches params.thisALG
        algrbhdf = algrbhdf[algrbhdf["gene_group"] == params.thisALG]
        #open the coo file
        coo = load_npz(input.coo)
        # In a previous version of this file, I did different things to replace the 0s with 9,999,999,999s in the missing data.
        # Now, I just leave everything as it is and only bother with the distances in the upstream script, PhyloTreeUMAP.snakefile

        # get the rows that have values that aren't in the index of algrbhdf
        # convert this to a pandas df with indices and colnames
        coodf = pd.DataFrame(coo.todense())
        # remove the rows that are not in the index of algrbhdf
        coodf = coodf.loc[algrbhdf.index]
        # remove the columns that are not in the index of algrbhdf
        coodf = coodf[algrbhdf.index]
        # change the indices to the value of the rbh column in the algrbhdf, matching on indices
        coodf.index   = [algrbhdf.loc[x, "rbh"] for x in coodf.index]
        # change the columns to the value of the rbh column in the algrbhdf, matching on columns
        coodf.columns = [algrbhdf.loc[x, "rbh"] for x in coodf.columns]
        #save the df to a tsv.gz file
        coodf.to_csv(output.filt_coo, sep="\t", compression="gzip", index=True)

def tsvgz_plotumap(sample, sampledffile, algrbhfile, tsvgz,
                   outdir, smalllargeNaN, n_neighbors, min_dist,
                   UMAPdf, UMAPbokeh):
    """
    This all-in-one plotting method makes UMAPs for the locus distance ALGs
        constructed by averaging across multiple species.
    Specifically, this is used for plotting the one-dot-one-locus UMAP plots.
    """
    # check that the types are correct
    if type(n_neighbors) not in [int, float]:
        raise ValueError(f"The n_neighbors {n_neighbors} is not of type int or float. Exiting.")
    if type(min_dist) not in [float]:
        raise ValueError(f"The min_dist {min_dist} is not of type float. Exiting.")
    # save the UMAP as a bokeh plot
    if not UMAPbokeh.endswith(".html"):
        raise ValueError(f"The output file {outhtml} does not end with '.html'. Exiting.")

    # read in the sample dataframe. We will need this later
    df = pd.read_csv(tsvgz, sep = "\t", index_col = 0)
    # read in the algrbh as a pandasdf
    algrbhdf = parse_rbh(algrbhfile)
    # only keep the rows that are in the df's indices that match the algrbhdf's rbh column'
    algrbhdf = algrbhdf[algrbhdf["rbh"].isin(df.index)]

    # read in the sampledf
    sampledf = pd.read_csv(sampledffile, sep = "\t", index_col = 0)
    # We need a unique set of files for each of these
    # In every case, we must produce a .df file and a .bokeh.html file
    # Sometimes, there are errors with making the UMAP, in that if the graph is not fully connected UMAP will fail.
    # For this reason, we should anticipate this error, and figure out what to do if it happens.
    print(f"    PLOTTING - UMAP with {smalllargeNaN} missing vals, with n_neighbors = {n_neighbors}, and min_dist = {min_dist}")
    reducer = umap.UMAP(n_neighbors = n_neighbors, min_dist = min_dist)
    disconnected = False
    with warnings.catch_warnings(record=True) as w:
        # Ignore UserWarnings temporarily
        warnings.filterwarnings("ignore", category=UserWarning)
        # Your code that might raise the warning
        mapper = reducer.fit(df.to_numpy())
        # Check if any warning was generated
        if w:
            for warning in w:
                if issubclass(warning.category, UserWarning):
                    disconnected = True
                    print("Got the warning that the graph is not fully connected. This happens mostly in the case of clades with highly conserved genomes:", warning.message)
                    # You can further process or log the warning here if needed
    # now we push left
    # save the UMAP as a bokeh plot
    if disconnected:
        plot_title = f"(Disconnected) Topo UMAP of {sample} with {smalllargeNaN} missing vals, n_neighbors = {n_neighbors}, min_dist = {min_dist}"
    else:
        plot_title = f"Topo UMAP of {sample} with {smalllargeNaN} missing vals, n_neighbors = {n_neighbors}, min_dist = {min_dist}"

    print("The embedding after fitting is:")
    print(mapper)
    print(mapper.embedding_)
    #              ┓    •
    # ┓┏┏┳┓┏┓┏┓  ┏┓┃┏┓╋╋┓┏┓┏┓
    # ┗┻┛┗┗┗┻┣┛  ┣┛┗┗┛┗┗┗┛┗┗┫
    #        ┛   ┛          ┛
    hover_data = pd.DataFrame({"rbh_ortholog": [-1] * len(df),
                               "gene_group":   [-1] * len(df),
                               "color":        [-1] * len(df)
                               })
    hover_data.index = algrbhdf["rbh"]
    hover_data["rbh_ortholog"] = hover_data.index
    # use apply to match the index to the "rbh" column of algrbhdf, and return the gene_group
    hover_data["gene_group"] = [algrbhdf[algrbhdf["rbh"] == x]["gene_group"].values[0] for x in hover_data.index]
    hover_data["color"]      = [algrbhdf[algrbhdf["rbh"] == x]["color"].values[0]      for x in hover_data.index]
    hover_data = hover_data.reset_index(drop = True)
    print(hover_data)
    # the index is now the name of the ortholog, so we have to get the colors another way
    #color_dict = dict(zip(algrbhdf["rbh"], algrbhdf["color"]))
    color_dict = dict(zip(hover_data.index, hover_data["color"]))
    print(color_dict)
    print(f"type of mapper is: {type(mapper)}")
    print(f"type of color_dict is: {type(color_dict)}")
    print(f"type of hover_data is: {type(hover_data)}")
    print(f"type of labels is: {type(list(hover_data['gene_group']))}")
    print(f"mapper is: \n", mapper)
    print(f"color_dict is: \n", color_dict)
    print(f"hover_data is: \n", hover_data)
    print(f"labels is: \n", list(hover_data["gene_group"]))

    # Needs to be umap 0.5.5 or later
    plot = umap.plot.interactive(mapper,
                                 color_key = color_dict,
                                 hover_data = hover_data,
                                 labels = list(hover_data["rbh_ortholog"]), # TODO this needs to be changd to a list comprehension
                                 tools=[], # this needs to be deleted, or else the zoom tool will not work.
                                 point_size = 4
                                 )
    # add a title to the plot
    plot.title.text = plot_title
    # output to an HTML file
    bokeh.io.output_file(UMAPbokeh)
    # Save the plot to an HTML file
    bokeh.io.save(plot)

    # get the coordinates of the UMAP
    df_embedding = pd.DataFrame(mapper.embedding_, columns=['UMAP1', 'UMAP2'])
    df_embedding.index = df.index
    umap_df = pd.concat([df, df_embedding], axis = 1)
    umap_df["gene_group"] = [algrbhdf[algrbhdf["rbh"] == x]["gene_group"].values[0] for x in umap_df.index]
    umap_df["color"]      = [algrbhdf[algrbhdf["rbh"] == x]["color"].values[0]      for x in umap_df.index]

    umap_df.to_csv(UMAPdf, sep = "\t", index = True)

def reimbedding_get_mem_mb(wildcards, attempt):
    """
    The amount of RAM needed could change.
    """
    attemptdict = {1: 2000,
                   2: 4000,
                   3: 8000,
                   4: 16000,
                   5: 32000,
                   6: 64000}
    return attemptdict[attempt]

rule ALG_reimbedding:
    input:
        tsvgz      = basedir + "/filt_coo/{analysis}_{ALG}.matrix.filt.tsv.gz",
        algrbhfile = config["ALG_rbh_file"],
        sampledf   = lambda wildcards: analysis_to_sampledf[wildcards.analysis]
    output:
        UMAPdf    = basedir + "/ALG_reimbedding/{analysis}/{ALG}/{analysis}.neighbors_{n}.mind_{m}.missing_large.subchrom.df",
        UMAPbokeh = basedir + "/ALG_reimbedding/{analysis}/{ALG}/{analysis}.neighbors_{n}.mind_{m}.missing_large.subchrom.bokeh.html"
    params:
        outdir   = lambda wildcards: basedir + "/ALG_reimbedding/{}/{}".format(wildcards.analysis, wildcards.ALG),
        analysis = lambda wildcards: wildcards.analysis,
        n        = lambda wildcards: int(wildcards.n),
        m        = lambda wildcards: float(wildcards.m)
    retries: 6
    threads: 1
    resources:
        mem_mb = reimbedding_get_mem_mb,
        runtime = 10
    run:
        tsvgz_plotumap(params.analysis, input.sampledf, input.algrbhfile, input.tsvgz,
                       params.outdir, "large", params.n, params.m,
                       output.UMAPdf, output.UMAPbokeh)

rule annotate_blastresults:
    """
    Annotates the UMAP df with the blast results, to make a plottable UMAP.
    Right now we don't do anything where we weight the point size phylogenetically.
    """
    input:
        UMAPdf     = basedir + "/ALG_reimbedding/{analysis}/{ALG}/{analysis}.neighbors_{n}.mind_{m}.missing_large.subchrom.df",
        algrbhfile = config["ALG_rbh_file"],
        blastp     = basedir + "/blast_concat/{query}/{analysis}.concat.filt.blastp",
    output:
        UMAPdf   = basedir + "/ALG_reimbedding/{analysis}/{ALG}/{analysis}.neighbors_{n}.mind_{m}.missing_large.subchrom.{query}.df"
    retries: 6
    threads: 1
    resources:
        mem_mb = tophit_get_mem_mb,
        runtime = 2
    run:
        # if the input UMAPdf is empty, write an empty df
        if os.stat(input.UMAPdf).st_size == 0:
            with open(output.UMAPdf, "w") as f:
                f.write("")
        else:
            # read in the df file
            UMAPdf = pd.read_csv(input.UMAPdf, sep = "\t", index_col = 0)
            UMAPdf["index"] = UMAPdf.index
            # read in the blastp results as a pandas df
            # If the file is empty, write an empty file to the output
            if os.stat(input.blastp).st_size == 0:
                with open(output.UMAPdf, "w") as f:
                    f.write("")
            else:
                # There are blast results, which likely means that the genome had an annotation and we didn't give it an Ersatz annotation with BCnS ALGs.
                # Keep going.
                blastp = pd.read_csv(input.blastp, sep = "\t", header = None)
                # make columns
                blastp.columns = ["sample", "qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore", "scaf_of_nearest_ALG", "nearest_ALG", "position_of_nearest_ALG", "dist_to_nearest_ALG"]
                blastp["closest_ALG"] = None
                blastp["closest_ALG_dist"] = None
                blastp["closest_ALG_position"] = None
                # get the blastp rows with nearest_ALG values in the UMAPdf.index
                for i, row in blastp.iterrows():
                    nearest_list = literal_eval(row["nearest_ALG"])
                    for j in range(len(nearest_list)):
                        if nearest_list[j] in UMAPdf.index:
                            nearest_pos = literal_eval(row["position_of_nearest_ALG"])[j]
                            nearest_dist = literal_eval(row["dist_to_nearest_ALG"])[j]
                            blastp.at[i, "closest_ALG"]          = nearest_list[j]
                            blastp.at[i, "closest_ALG_dist"]     = nearest_dist
                            blastp.at[i, "closest_ALG_position"] = nearest_pos
                            break
                # remove rows where closest_ALG is None
                blastp = blastp[blastp["closest_ALG"].notnull()]

                # groupby the qseqid column
                gb = blastp.groupby("qseqid")
                newcols = {}
                for name, group in gb:
                    # add the number of times that this index appears in this group
                    newcols[f"PLOTGENE_{name}"] = UMAPdf["index"].apply(lambda x: group[group["closest_ALG"] == x].shape[0])
                # add the new columns to the UMAPdf
                newcolsdf = pd.DataFrame(newcols)
                UMAPdf = pd.concat([UMAPdf, newcolsdf], axis = 1)
                # sum up all the counts across all query genes
                UMAPdf["plotgene_SUM"] = UMAPdf[[x for x in UMAPdf.columns if "PLOTGENE_" in x]].sum(axis = 1)

                algrbhdf = parse_rbh(input.algrbhfile)
                UMAPdf["color"] = [algrbhdf[algrbhdf["rbh"] == x]["color"].values[0] for x in UMAPdf.index]

                # remove the "index" column
                UMAPdf = UMAPdf.drop("index", axis = 1)
                print(list(UMAPdf.columns))
                # add a color column
                UMAPdf.to_csv(output.UMAPdf, sep = "\t", index = True)
                # A joke for you
                # Why did the chicken cross the road?
                # To get to the other side where the programmer was coding a virtual chicken-crossing algorithm in Python!
                #   - Attribution: Chat-GPT v3.5, Monday March 11 5:55PM

def umapdf_to_pdf(UMAPdf, analysis, ALG, n, m, query, outputPDF, species = None):
    """
    This function takes a dataframe as input and makes a UMAP
    """
    #         ┓   ┓•┓
    # ┏┳┓┏┓╋┏┓┃┏┓╋┃┓┣┓
    # ┛┗┗┗┻┗┣┛┗┗┛┗┗┗┗┛
    #       ┛
    # If the UMAPdf is empty, write a plot telling the user that the input df was empty, so we couldn't plot anything
    if os.stat(UMAPdf).st_size == 0:
        # make a 2in x 2in plot telling the user the message
        fig = plt.figure(figsize=(2,2))
        plt.text(0.5, 0.5, "The input df was empty, so we couldn't plot anything", fontsize = 3)
        # make another line telling the user which file, exactly, was empty
        plt.text(0.5, 0.4, f"The input file was {UMAPdf}", fontsize = 3)
        # turn off the axis ticks
        for ax in fig.axes:
            ax.set_xticks([])
            ax.set_yticks([])
        # turn off the axes
        plt.axis('off')
        try:
            # make the plot tight to not cut off the text
            plt.tight_layout()
        except:
            pass
        # save the figure
        plt.savefig(outputPDF)
    else:
        dot = 3
        bigger_dot = 4
        df_embedding = pd.read_csv(UMAPdf, sep = "\t", index_col = 0)
        # make a matplotlib plot of the UMAP with the df_embedding, and the color_dict from SplitLossColocTree as the legend
        # make a figure that is 5x5 inches
        fig = plt.figure(figsize=(2,2))
        # scatter the UMAP1 and UMAP2 columns of the df_embedding
        plt.scatter(df_embedding["UMAP1"], df_embedding["UMAP2"], c = df_embedding["color"], s = dot)
        subdf = df_embedding[df_embedding["plotgene_SUM"] > 0]
        plt.scatter(subdf["UMAP1"], subdf["UMAP2"], c = "#5AFF00", s = [dot + (dot*bigger_dot*(x/subdf["plotgene_SUM"].max())) for x in subdf["plotgene_SUM"]], alpha = 0.25)
        ## make a legend with the color_dict from SplitLossColocTree
        #nbci = NCBITaxa()
        ## get the name of the ncbi taxid from the SplitLossColocTree color_dict
        #legend_dict = {}
        #for key in SplitLossColocTree.color_dict_top:
        #    taxid = int(key)
        #    taxname = nbci.get_taxid_translator([taxid])[taxid]
        #    legend_dict[taxname] = SplitLossColocTree.color_dict_top[key]
        #print("This is the legend dict")
        #print(legend_dict)
        #legend_patches = [mpatches.Patch(color=color, label=label)
        #                  for label, color in legend_dict.items()]
        ## add the entries to the legend
        #fig = plt.legend(handles=legend_patches, loc="upper right")
        # turn off the axis ticks
        for ax in fig.axes:
            ax.set_xticks([])
            ax.set_yticks([])
        # set the title based on the input
        plt.title(f"UMAP of {analysis}, ALG {ALG},\nwith {query} blast results,\nmin_dist {m}, n_neighbors {n}",
                   fontsize = 3)
        # save the figure
        plt.savefig(outputPDF)

rule plot_clade_blast:
    input:
        UMAPdf   = basedir + "/ALG_reimbedding/{analysis}/{ALG}/{analysis}.neighbors_{n}.mind_{m}.missing_large.subchrom.{query}.df",
    output:
        outPDF   = basedir + "/ALG_reimbedding/{analysis}/{ALG}/{analysis}.neighbors_{n}.mind_{m}.missing_large.subchrom.{query}.pdf"
    retries: 6
    threads: 1
    resources:
        mem_mb = tophit_get_mem_mb,
        runtime = 1
    run:
        umapdf_to_pdf(input.UMAPdf, wildcards.analysis, wildcards.ALG, wildcards.n, wildcards.m, wildcards.query, output.outPDF)


rule plot_one_species_one_query:
    """
    The goal of this rule is to plot the UMAP of a specific ALG for a specific query.
    The UMAP embedding is from the analysis to which the species belongs.
    Within this plot, we assign a color to each ALG, and also make a legend to see which genes are where.
    """
    input:
        UMAPdf = basedir + "/ALG_reimbedding/{analysis}/{ALG}/{analysis}.neighbors_{n}.mind_{m}.missing_large.subchrom.df",
        blastp = basedir + "/blast_filt/{query}/{sample}_results.filt.blastp",
        algrbhfile = config["ALG_rbh_file"],
    output:
        outPDF = basedir + "/ALG_reimbedding/one_analysis_one_species_one_query/{analysis}.{sample}.{query}.{ALG}.neighbors_{n}.mind_{m}.missing_large.subchrom.pdf"
    threads: 1
    resources:
        mem_mb = tophit_get_mem_mb,
        runtime = 1
    run:
        umapdf_one_species_one_query(input.UMAPdf, input.blastp, wildcards.analysis,
                                     wildcards.ALG, wildcards.n, wildcards.m,
                                     wildcards.query, output.outPDF, species = wildcards.sample)

#def plot_paramsweep_one_species_one_query(dataframe_files, blastp):
#    """
#    Makes the plot for the parameter sweep plot when we provide multiple dataframes.
#    """
#    df_filelist = []
#    if args.directory:
#        # Get all the dataframes from the directory
#        df_filelist = [x for x in os.listdir(args.directory)
#                       if x.endswith(".df")]
#    elif args.filelist:
#        # Get all the dataframes from the filelist
#        df_filelist = args.filelist.split(" ")
#
#    print("df_filelist is {}".format(df_filelist))
#
#    # There filenames are formatted like this:
#    #   Metazoa_33208.neighbors_20.mind_0.0.missing_large.df
#    samplename = set([os.path.basename(x).split(".neighbors_")[0] for x in df_filelist])
#    if len(samplename) > 1:
#        raise ValueError("More than one sample name found in the file list. We only allow one for now.")
#    samplename = list(samplename)[0]
#
#    print("samplename is {}".format(samplename))
#
#    # get the number of neighbors
#    num_neighbors = [os.path.basename(x).split(".neighbors_")[1].split(".")[0] for x in df_filelist]
#    # check that everything can be cast to an int
#    for x in num_neighbors:
#        try:
#            int(x)
#        except ValueError:
#            raise ValueError(f"Could not cast {x} to an int.")
#    num_neighbors = [int(x) for x in num_neighbors]
#    print("num_neighbors is {}".format(num_neighbors))
#    # get the min_dist
#    min_dist = [os.path.basename(x).split(".mind_")[1].split(".missing_")[0] for x in df_filelist]
#    # check that all min_dist can be cast to a floar
#    for x in min_dist:
#        try:
#            float(x)
#        except ValueError:
#            raise ValueError(f"Could not cast {x} to a float.")
#    min_dist = [float(x) for x in min_dist]
#    print("min_dist is {}".format(min_dist))
#    # get whether it is from the small or large dataset
#    miss_size = [os.path.basename(x).split(".missing_")[1].split(".")[0] for x in df_filelist]
#    # make sure that there is only one value for small or large here. We can't deal with both.
#    if len(set(miss_size)) > 1:
#        raise ValueError("More than one size found in the file list. We only allow one for now: small or large")
#    print("Missing size is {}".format(miss_size))
#
#    # collate together the html_filelist, num_neighbors, min_dist, and size
#    df = pd.DataFrame({"dffile":        df_filelist,
#                       "num_neighbors": num_neighbors,
#                       "min_dist":      min_dist,
#                       "size":          miss_size})
#    # sort by num_neighbors and min_dist, ascending both
#    df = df.sort_values(["num_neighbors", "min_dist"], ascending=[True, True])
#    # groupby size and make one plot for each size
#    gb = df.groupby("size")
#    for name, group in gb:
#        # make subplots.
#        # The number of unique things in num_neighbors is the number of rows
#        # the numer of unique things in min_dist is the number of columns
#        num_rows = len(group["num_neighbors"].unique())
#        num_cols = len(group["min_dist"].unique())
#        # make a grid of squares to plot each of these on
#        # reduce the space between all the plots
#        # make the figure size such that all the plots are square
#        square_size = 1.5
#        fig, axes = plt.subplots(num_rows, num_cols,
#                                 figsize=(num_cols*square_size, num_rows*square_size))
#        plt.subplots_adjust(wspace=0.1, hspace=0.1)
#
#        # turn off all the ticks
#        # Turn off the axes
#        for i in range(num_rows):
#            for j in range(num_cols):
#                axes[i, j].set_xticks([])
#                axes[i, j].set_yticks([])
#                # Turn off the lines around the plot
#                for spine in ['top', 'right', 'bottom', 'left']:
#                    axes[i, j].spines[spine].set_visible(False)
#
#        # set the title as samplename and the whether it is small or large
#        fig.suptitle(f"{samplename} {name} NaNs")
#        # set absolute left label as the number of neighbors
#        fig.text(0.06, 0.5, 'Number of Neighbors', va='center', rotation='vertical')
#        fig.text(0.5, 0.92, 'Min Distance', ha='center')
#        for i, row in enumerate(sorted(group["num_neighbors"].unique(), reverse = False)):
#            # for the left-most plot in each row, set the ylabel as the number of neighbors
#            axes[i, 0].set_ylabel(f"{row}", rotation=0, ha='right')
#            for j, col in enumerate(sorted(group["min_dist"].unique(), reverse = False)):
#                #for the top-most plot in each column, set the xlabel as the min_dist
#                # put the xlabel on the top
#                axes[0, j].xaxis.set_label_position('top')
#                axes[0, j].set_xlabel(f"{col}")
#                # get the df file for this row and column from the dffile
#                dffile = group[(group["num_neighbors"] == row) & (group["min_dist"] == col)]["dffile"].values[0]
#                # if the type of dffile is NoneType, then we didn't find the file
#                if type(dffile) == type(None):
#                    # write into the ax[i, j] that we didn't find the file
#                    axes[i, j].text(0.5, 0.5, f"Missing file", fontsize=6, ha='center')
#                elif type(dffile) == str:
#                    # now we try to find the file
#                    if os.path.exists(dffile):
#                        # if the file is empty, then we write into the ax[i, j] that the file is empty
#                        if os.path.getsize(dffile) == 0:
#                            axes[i, j].text(0.5, 0.5, f"Empty file", fontsize=6, ha='center')
#                        else:
#                            tempdf = pd.read_csv(dffile, sep="\t", index_col=0)
#                            # Plot the UMAP1 UMAP2 with the color column as the color.
#                            # Make the dot size small
#                            axes[i, j].scatter(tempdf["UMAP1"], tempdf["UMAP2"],
#                                             s=0.5, lw = 0, alpha=0.5,
#                                             color=list(tempdf["color"]))
#                    else:
#                        # The user provided the path to this file, but it doesn't exist.
#                        # This means that the user made a mistake in writing the file name.
#                        raise ValueError(f"The file {dffile} does not exist.")
#                else:
#                    raise ValueError(f"Type of dffile is not a string or NoneType. It is {type(dffile)}")
#        # make sure that the aspect ratio for all of these is the same
#        # Iterate over each axis and set aspect ratio to 'equal'
#        for row in axes:
#            for col in row:
#                col.set_aspect('equal', adjustable='box')
#
#        # Now make vertical and horizontal lines to separate the plots. Make them medium gray.
#        # The lines will be on the figure, and not in the plots.
#        # Draw horizontal lines between rows
#        # Get the bounding boxes of the axes including text decorations
#
#        # Get the bounding boxes of the axes including text decorations
#        r = fig.canvas.get_renderer()
#        get_bbox = lambda ax: ax.get_tightbbox(r).transformed(fig.transFigure.inverted())
#        bbox_list = [get_bbox(ax) for ax in axes.flat]
#
#        # Create an empty array with the correct shape and dtype
#        bboxes = np.empty(axes.shape, dtype=object)
#
#        # Fill the array with the bounding boxes
#        for idx, bbox in np.ndenumerate(bboxes):
#            bboxes[idx] = bbox_list[idx[0] * axes.shape[1] + idx[1]]
#
#        # Get the minimum and maximum extent, get the coordinate half-way between those
#        ymax = np.array(list(map(lambda b: b.y1, bboxes.flat))).reshape(axes.shape).max(axis=1)
#        ymin = np.array(list(map(lambda b: b.y0, bboxes.flat))).reshape(axes.shape).min(axis=1)
#        ys = np.c_[ymax[1:], ymin[:-1]].mean(axis=1)
#
#        # Draw horizontal lines at those coordinates
#        for y in ys:
#            line = plt.Line2D([0.125, 0.9], [y, y], transform=fig.transFigure, color="#BBBBBB")
#            fig.add_artist(line)
#
#        # Get the minimum and maximum extent, get the coordinate half-way between those for vertical lines
#        xmax = np.array(list(map(lambda b: b.x1, bboxes.flat))).reshape(axes.shape).max(axis=0)
#        xmin = np.array(list(map(lambda b: b.x0, bboxes.flat))).reshape(axes.shape).min(axis=0)
#        xs = np.c_[xmax[1:], xmin[:-1]].mean(axis=1)
#
#        # Draw vertical lines at those coordinates
#        for xi in range(len(xs)):
#            x = xs[xi]
#            if xi == 0:
#                x = x + 0.0125
#            line = plt.Line2D([x, x], [0.1, 0.875], transform=fig.transFigure, color="#BBBBBB")
#            fig.add_artist(line)
#
#        # save the figure as f"{samplename}_{name}.pdf"
#        # name is just the size, small or large
#        print("saving the file to {}".format(args.outpdf))
#        plt.savefig(args.outpdf)
#        # close the figure
#        plt.close(fig)
#
#rule paramsweep_one_species_one_query:
#    """
#    This rule is for making a plot of the different parameters that are plotted for a single species,
#    and a single query.
#    """
#    input:
#        UMAPdf = expand(basedir + "/ALG_reimbedding/{{analysis}}/{{ALG}}/{{analysis}}.neighbors_{n}.mind_{m}.missing_large.subchrom.df",
#                n = odol_n,
#                m = odol_m),
#        blastp = basedir + "/blast_filt/{query}/{sample}_results.filt.blastp",
#    output:
#        sweep = basedir + "/ALG_reimbedding/one_analysis_one_species_one_query/sweeps/{{analysis}}.{{sample}}.{{query}}.{{ALG}}.missing_large.subchrom.sweep.pdf"
#    output:
