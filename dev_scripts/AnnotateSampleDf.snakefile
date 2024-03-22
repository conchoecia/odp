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

ofix = "dfannotate"

rule all:
    input:
        ofix + "/measurements/allsamples.protstats.collated.df",
        ofix + "/measurements/allsamples.rbhstats.collated.df",
        ofix + "/measurements/allsamples.genomestats.collated.df",
        expand(ofix + "/{dfname}.supplemented.df",
            dfname=dfname_to_filepath.keys()),
        expand(ofix + "/{dfname}.supplemented.pdf",
            dfname=dfname_to_filepath.keys())

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