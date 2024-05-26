"""
Program  : AnnotateSampleDf.snakefile
Language : snakemake
Date     : 2024-02-29
Author   : Darrin T. Schultz
Email    : darrin.schultz@univie.ac.at
Github   : https://github.com/conchoecia/odp
Support  : For issues or questions, please search if the topic has been discussed already
           on github and open a new issue if not: https://github.com/conchoecia/odp/issues
License  : GNU GENERAL PUBLIC LICENSE, Version 3, 29 June 2007. See the LICENSE file.
Citation : None currently

Description:
  This program reads in all the genomes from a config file and makes a dataframe of the stats.

  Also takes in a directory of rbh files and makes measurements from that.

  Takes in a rbh df to infer whether genomes are annotated by the RBH, or if they are from another source.

  Takes a list of input .df files and annotates them with the genome information.


Usage instructions:
  - See https://github.com/conchoecia/odp#getting-started
"""

# This block imports fasta-parser as fasta
import os
import sys
snakefile_path = os.path.dirname(os.path.realpath(workflow.snakefile))
dependencies_path = os.path.join(snakefile_path, "../dependencies/fasta-parser")
sys.path.insert(1, dependencies_path)
import fasta
# other imports
import pandas as pd
import AnnotateSampleDf as asd

configfile: "config.yaml"

# check that ALG_name is in the config
if not "ALG_name" in config:
    raise ValueError("ALG_name not in config")

# for every species in the config - make sure there is a RBH file
config["species_to_rbh"] = {}
for thisrbh in os.listdir(config["rbh_dir"]):
    if thisrbh.endswith(".rbh") and thisrbh.startswith(config["ALG_name"]):
        sample = thisrbh.split("_")[1]
        config["species_to_rbh"][sample] = os.path.join(config["rbh_dir"], thisrbh)

in_sample_not_in_rbh = set(config["species"].keys()) - set(config["species_to_rbh"].keys())
if len(in_sample_not_in_rbh) > 0:
    raise ValueError(f"Samples {in_sample_not_in_rbh} in species_to_rbh but not in species")

# make sure there is a dfs entry in config
if not "dfs" in config:
    raise ValueError("dfs not in config")
dfname_to_filepath = {os.path.basename(x).replace(".df","").replace(".tsv",""): x
                      for x in config["dfs"]}
print(dfname_to_filepath)


# If we want to make subplots of specific clades, execute the code below.
# There are a lot of checks to perform.
if "subplots" in config:
    for subplot in config["subplots"]:
        # make sure that each key in the subplots has at least "fileprefix", "title", and "taxids" keys
        for key in ["fileprefix", "title", "taxids_to_include"]:
            if not key in subplot:
                raise ValueError(f"Subplot {subplot} does not have a {key} key")
        # assert types
        if type(subplot["fileprefix"]) != str:
            raise ValueError(f"Subplot key for {subplot} is not a string")
        if type(subplot["title"]) != str:
            raise ValueError(f"Subplot key for {subplot} is not a string")
        if type(subplot["taxids_to_include"]) != list:
            raise ValueError(f"Subplot key for {subplot} is not a list")
        # check that the fileprefixes of the subplots are all just alphanumeric
        if not subplot["fileprefix"].isalnum():
            raise ValueError(f"Subplot key for {subplot} is not alphanumeric. We only allow the chars [a-zA-Z0-9]")
        # Check that the values in "taxids_to_include" are all integers
        for thisval in subplot["taxids_to_include"]:
            if not isinstance(thisval, int):
                raise ValueError(f"Value {thisval} in subplot key {subplot} is not an integer. We only allow integers")
        # Check all of the keys, and only allow the three we mentioned, plus "taxids_to_exclude"
        for thiskey in subplot.keys():
            if thiskey not in ["fileprefix", "title", "taxids_to_include", "taxids_to_exclude"]:
                raise ValueError(f"Subplot key for {subplot} is not allowed. We only allow fileprefix, title, taxids, and taxids_to_include")

        # if taxids_to_exclude is in the dictionary, make sure that all of the values are ints
        if "taxids_to_exclude" in subplot:
            if type(subplot["taxids_to_exclude"]) != list:
                raise ValueError(f"Subplot key for {subplot} is not a list")
            for thisval in subplot["taxids_to_exclude"]:
                if not isinstance(thisval, int):
                    raise ValueError(f"Value {thisval} in {subplot} is not an integer. We only allow integers")
        # now that we are sure that subplots is legal, reformat such that the fileprefix is the key to one big dict
    config["subplots"] = {s["fileprefix"]: s for s in config["subplots"]}

print(config["subplots"])

ofix = "dfannotate"

rule all:
    input:
        #ofix + "/measurements/allsamples.protstats.collated.df",
        #ofix + "/measurements/allsamples.rbhstats.collated.df",
        #ofix + "/measurements/allsamples.genomestats.collated.df",
        #expand(ofix + "/{dfname}.supplemented.df",
        #    dfname=dfname_to_filepath.keys()),
        #expand(ofix + "/{dfname}.supplemented.pdf",
        #    dfname=dfname_to_filepath.keys()),
        ## this is the ALG dispersion plot
        #ofix + "/measurements/rbh_dispersion_plot.pdf",
        ## these are all the subclade plots
        #expand(ofix + "/subplots/{dfname}_{subplot}.pdf",
        #       dfname = dfname_to_filepath.keys(),
        #       subplot = config["subplots"].keys()),
        # also produce the data availability statement for all the genomes
        expand(ofix + "/dataAvailability_{dfname}.html",
               dfname = dfname_to_filepath.keys())

rule data_availability:
    """
    This rule reads in the dataframe and produces a Data Availability statement, complete with URLs.
    """
    input:
        df = lambda wildcards: dfname_to_filepath[wildcards.dfname]
    output:
        html = ofix + "/dataAvailability_{dfname}.html"
    threads: 1
    retries: 1
    resources:
        mem_mb  = 1000,
        runtime = 5
    run:
        asd.gen_data_availability_statement(input.df, output.html)


def stats_get_mem_mb(wildcards, attempt):
    """
    The amount of RAM needed for the script depends on the size of the input genome.
    """
    attemptdict = {1: 250,
                   2: 1000,
                   3: 5000,
                   4: 20000,
                   5: 40000,
                  }
    return attemptdict[attempt]

def stats_get_runtime(wildcards, attempt):
    """
    The amount of RAM needed for the script depends on the size of the input genome.
    """
    attemptdict = {1: 2,
                   2: 5,
                   3: 10,
                   4: 20,
                   5: 40,
                  }
    return attemptdict[attempt]

rule genstats:
    """
    This generates statistics about the genome, and saves them to key: value pairs in a text file.
    The fields that are output are:
      - the number of scaffolds
      - the GC content
      - the genome size
      - the median scaffold length
      - the mean scaffold length
      - scaffold N50
      - longest scaffold
      - smallest scaffold
      - percent Ns
    """
    input:
        genome = lambda wildcards: config["species"][wildcards.sample]["genome"]
    output:
        results = ofix + "/measurements/genome_stats/{sample}.genome_stats.txt"
    threads: 1
    retries: 5
    resources:
        mem_mb = stats_get_mem_mb,
        runtime   = stats_get_runtime
    run:
        asd.gen_genome_stats(input.genome, output.results)

rule annotation_stats:
    """
    This generates statistics about the protein composition of the genome.
    The information that is saved is:
      - the number of proteins
      - the mean protein length
      - the median protein length
      - the longest protein
      - the smallest protein
      - whether the proteins are from a real annotation or from the RBH entries
    """
    input:
        proteins = lambda wildcards: config["species"][wildcards.sample]["proteins"],
        alg_rbh  = config["ALG_rbh"]
    output:
        results = ofix + "/measurements/protein_stats/{sample}.protein_stats.txt"
    threads: 1
    retries: 5
    resources:
        mem_mb = stats_get_mem_mb,
        runtime   = stats_get_runtime
    run:
        asd.gen_annotation_stats(input.proteins, input.alg_rbh, output.results)

rule rbhstats:
    """
    This rule generates statistics about the rbh files.
    Namely, it characterizes the dispersion properties.
    - The fields that are output are:
      - frac_ologs:           The fraction of genes of ANY ALG that are present at all in the rbh file. len(rbhdf) / total_genes_ALGs
      - frac_ologs_sig:       The fraction of genes of ANY ALG that are significantly on any chromosome, as defined by whole_FET
      - frac_ologs_single:    The fraction of genes of ANY ALG that are significantly on the largest chromosome, as defined by whole_FET
      - frac_ologs_{ALGNAME}: The fraction of genes of INDIVIDUAL ALGs that are significantly on any chromosome
    """
    input:
        rbh = lambda wildcards: config["species_to_rbh"][wildcards.sample],
        alg_rbh  = config["ALG_rbh"]
    output:
        results = ofix + "/measurements/rbh_stats/{sample}.rbh_stats.txt"
    threads: 1
    retries: 5
    resources:
        mem_mb  = stats_get_mem_mb,
        runtime = stats_get_runtime
    run:
        asd.gen_rbh_stats(input.rbh, input.alg_rbh, config["ALG_name"], output.results)

def composite_get_mem_mb(wildcards, attempt):
    """
    The amount of RAM needed for the script depends on the size of the input genome.
    """
    attemptdict = {1: 1000,
                   2: 5000,
                   3: 10000,
                   4: 20000,
                   5: 40000,
                  }
    return attemptdict[attempt]

def composite_get_runtime(wildcards, attempt):
    """
    The amount of RAM needed for the script depends on the size of the input genome.
    """
    attemptdict = {1: 10,
                   2: 20,
                   3: 40,
                   4: 80,
                   5: 160,
                  }
    return attemptdict[attempt]

rule CompProt:
    """
    make a composite file for the proteins
    """
    input:
        prot_stats = expand(ofix + "/measurements/protein_stats/{sample}.protein_stats.txt",
                            sample=config["species"]),
    output:
        df = ofix + "/measurements/allsamples.protstats.collated.df"
    threads: 1
    retries: 5
    resources:
        mem_mb  = composite_get_mem_mb,
        runtime = composite_get_runtime
    run:
        # In this rule, we merge together the dataframe with the genome stats, the protein stats, and the rbh stats
        entries = []
        for sample in config["species"]:
            mydict = {"sample": sample}
            # read in the prot_stats
            prot_stats_file   = f"{ofix}/measurements/protein_stats/{sample}.protein_stats.txt"
            mydict.update(asd.stats_filepath_to_dict(prot_stats_file))
            entries.append(mydict)
        updatedf = pd.DataFrame(entries)
        # save to the output file
        updatedf.to_csv(output.df, sep="\t")

rule CompRBH:
    """
    make a composite file for the info from the rbh files
    """
    input:
        rbh_stats = expand(ofix + "/measurements/rbh_stats/{sample}.rbh_stats.txt",
                            sample=config["species"])
    output:
        df = ofix + "/measurements/allsamples.rbhstats.collated.df"
    threads: 1
    retries: 5
    resources:
        mem_mb  = composite_get_mem_mb,
        runtime = composite_get_runtime
    run:
        # In this rule, we merge together the dataframe with the genome stats, the protein stats, and the rbh stats
        entries = []
        for sample in config["species"]:
            mydict = {"sample": sample}
            # read in the prot_stats
            rbh_stats_file    = f"{ofix}/measurements/rbh_stats/{sample}.rbh_stats.txt"
            mydict.update(asd.stats_filepath_to_dict(rbh_stats_file))
            entries.append(mydict)
        updatedf = pd.DataFrame(entries)
        # save to the output file
        updatedf.to_csv(output.df, sep="\t")

rule CompGenome:
    """
    make a composite file for the info from the genome fasta files
    """
    input:
        genome_stats = expand(ofix + "/measurements/genome_stats/{sample}.genome_stats.txt",
                            sample=config["species"]),
    output:
        df = ofix + "/measurements/allsamples.genomestats.collated.df"
    threads: 1
    retries: 5
    resources:
        mem_mb  = composite_get_mem_mb,
        runtime = composite_get_runtime
    run:
        # In this rule, we merge together the dataframe with the genome stats, the protein stats, and the rbh stats
        entries = []
        for sample in config["species"]:
            mydict = {"sample": sample}
            # read in the prot_stats
            genome_stats_file = f"{ofix}/measurements/genome_stats/{sample}.genome_stats.txt"
            mydict.update(asd.stats_filepath_to_dict(genome_stats_file))
            entries.append(mydict)
        updatedf = pd.DataFrame(entries)
        # save to the output file
        updatedf.to_csv(output.df, sep="\t")

rule make_composite_dataframe:
    input:
        df       = lambda wildcards: dfname_to_filepath[wildcards.dfname],
        g_stats = ofix + "/measurements/allsamples.genomestats.collated.df",
        r_stats = ofix + "/measurements/allsamples.rbhstats.collated.df",
        p_stats = ofix + "/measurements/allsamples.protstats.collated.df",
    output:
        df = ofix + "/{dfname}.supplemented.df"
    threads: 1
    resources:
        mem_mb  = 250,
        runtime = 2
    run:
        # import the three dataframes to join by sample column
        compositedf = pd.read_csv(input.df, sep="\t", index_col=0)
        for df in [input.r_stats, input.g_stats, input.p_stats]:
            # do a left merge on the sample column in both dfs
            tempdf = pd.read_csv(df, sep="\t", index_col=0)
            # if there are columns in compositedf that already exist in tempdf, drop them from compositedf
            compositedf = compositedf.drop(columns=[col for col in compositedf.columns if (col in tempdf.columns) and (col != "sample")])
            # merge on the sample column
            compositedf = compositedf.merge(tempdf, left_on="sample", right_on="sample", how="left")
        compositedf.to_csv(output.df, sep="\t", index=True)

rule pdf:
    """
    This makes a plot of all the statistics, with each statistic plotted over the samples
    """
    input:
        df = ofix + "/{dfname}.supplemented.df",
        plotdfs = os.path.join(snakefile_path, "PhyloTreeUMAP_plotdfs.py")
    output:
        pdf = ofix + "/{dfname}.supplemented.pdf"
    threads: 1
    resources:
        mem_mb  = 1000,
        runtime = 3
    shell:
        """
        python {input.plotdfs} --plot_features -f {input.df} -o {output.pdf}
        """

rule dispersion_plot:
    """
    This generates a dispersion plot of the rbh files from a single dataframe.
    """
    input:
        df = ofix + "/measurements/allsamples.rbhstats.collated.df",
        alg_rbh  = config["ALG_rbh"]
    output:
        pdf = ofix + "/measurements/rbh_dispersion_plot.pdf"
    threads: 1
    resources:
        mem_mb  = 1000,
        runtime = 3
    run:
        asd.bin_and_plot_decay(input.alg_rbh, input.df, output.pdf, config["ALG_name"], 5)

rule subclade_plots:
    """
    These are plots of the subclades that are defined in the config file.
    The output of this file is a pdf that has the dots of the specified clade highlighted, and everything else is not highlighted.
    """
    input:
        df  = ofix + "/{dfname}.supplemented.df"
    output:
        pdf = ofix + "/subplots/{dfname}_{subplot}.pdf"
    threads: 1
    resources:
        mem_mb  = 1000,
        runtime = 5
    params:
        title = lambda wildcards: config["subplots"][wildcards.subplot]["title"],
        taxids_to_include = lambda wildcards: config["subplots"][wildcards.subplot]["taxids_to_include"],
        taxids_to_exclude = lambda wildcards: config["subplots"][wildcards.subplot].get("taxids_to_exclude", [])
    run:
        asd.plot_UMAP_highlight_subclade(input.df,
                                         params.title,
                                         params.taxids_to_include,
                                         params.taxids_to_exclude,
                                         output.pdf)
