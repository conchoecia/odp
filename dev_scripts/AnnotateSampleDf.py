#!/usr/bin/env python

"""
This program contains the functions used by AnnotateSampleDf.snakefile
"""

# This block imports fasta-parser as fasta
import os
import sys
thispath = os.path.dirname(os.path.realpath(__file__))
dependencies_path = os.path.join(thispath, "../dependencies/fasta-parser")
sys.path.insert(1, dependencies_path)
import fasta

import pandas as pd
from rbh_tools import parse_rbh

def gen_rbh_stats(samplerbhfilepath, algrbhfilepath, ALGname, outfilepath):
    """
    This function generates the stats of an rbh file - namely the dispersion.
    Things that are calculated for this are:
     - the number of proteins in the rbh file
     - the number of gene groups in the rbh file
     - the number of genes for each gene group

    Input:
      - It takes in a single argument, the path to one rbh file.
    Output:
      - The output is a text file that contains the analysis as key: value pairs.
      - The fields that are output are:
        - frac_ologs:           The fraction of genes of ANY ALG that are present at all in the rbh file. len(rbhdf) / total_genes_ALGs
        - frac_ologs_sig:       The fraction of genes of ANY ALG that are significantly on any chromosome, as defined by whole_FET
        - frac_ologs_single:    The fraction of genes of ANY ALG that are significantly on the largest chromosome, as defined by whole_FET
        - frac_ologs_{ALGNAME}: The fraction of genes of INDIVIDUAL ALGs that are significantly on any chromosome
    """
    ALGrbhdf = parse_rbh(algrbhfilepath)
    total_genes_ALGs = len(ALGrbhdf)
    genes_per_ALG    = ALGrbhdf.groupby("gene_group").size().to_dict()

    # now parse the rbh file
    rbhdf = parse_rbh(samplerbhfilepath)
    sigdf = rbhdf[rbhdf["whole_FET"] < 0.05]

    # get the sample scaf column
    sample_loc_col = [col for col in rbhdf.columns if (col.endswith("_scaf")) and (ALGname not in col)][0]
    frac_ologs = len(rbhdf)/total_genes_ALGs
    frac_ologs_sig = len(rbhdf[rbhdf["whole_FET"] < 0.05])/total_genes_ALGs
    # groupby the gene_group, then get the rows with the most frequent sample_loc_col value
    entries = []
    for gene_group, groupdf in sigdf.groupby("gene_group"):
        max_single = 0 if len(groupdf) == 0 else groupdf[sample_loc_col].value_counts().max()
        entries.append({"gene_group":     gene_group,
                        "max_single":     max_single,
                        "genes_in_group": genes_per_ALG[gene_group]})
    # frac_ologs_single is the sum o
    if len(entries) == 0:
        frac_ologs_single = float(0)
    else:
        frac_ologs_single = pd.DataFrame(entries)["max_single"].sum() / total_genes_ALGs
    # print all of these fields to a text file in the format:
    # field: value
    with open(outfilepath, "w") as f:
        f.write(f"frac_ologs: {frac_ologs}\n")
        f.write(f"frac_ologs_sig: {frac_ologs_sig}\n")
        f.write(f"frac_ologs_single: {frac_ologs_single}\n")
        for ALG in genes_per_ALG:
            f.write(f"frac_ologs_{ALG}: {len(sigdf[sigdf['gene_group'] == ALG])/genes_per_ALG[ALG]}\n")

def gen_annotation_stats(sampleproteinfilepath, algrbhfilepath, outfilepath):
    """
    This generates information about the genome annotation. Specifically, it looks at the proteins in the annotation.
    Things that are calculated:
      - the number of proteins
      - the mean protein length
      - the median protein length
      - the longest protein
      - the smallest protein
      - whether the proteins are from a real annotation or from the RBH entries
    Input:
      - algrbhfile
      - protein fasta.gz file
    output:
      - A text file with the following fields:
        - num_proteins: {num_proteins}
        - mean_protein_length: {mean_protein_length}
        - median_protein_length: {median_protein_length}
        - longest_protein: {longest_protein}
        - smallest_protein: {smallest_protein}
        - from_rbh: {from_rbh}
    """
    # read in the ALG_rbh file as a pandas df
    df = parse_rbh(algrbhfilepath)
    rbh_names = list(df["rbh"])
    # read in the proteins. Make a list of putative rbh proteins. Get the other stats.
    entries = []
    for record in fasta.parse(sampleproteinfilepath):
        entries.append({"protname":           record.id,
                        "protlen" :           len(record.seq),
                        "putative_rbh_name" : "_".join(record.id.split("_")[:-1]) if "_" in record.id else record.id })
    protdf = pd.DataFrame(entries)
    num_proteins = len(protdf)
    mean_protein_length = protdf["protlen"].mean()
    median_protein_length = protdf["protlen"].median()
    longest_protein  = protdf["protlen"].max()
    smallest_protein = protdf["protlen"].min()
    # count the number of times the putative_rbh_name is in the rbh_names
    if protdf["putative_rbh_name"].isin(rbh_names).sum() > (0.25 * len(protdf)):
        from_rbh = True
    else:
        from_rbh = False
    # Print all of these fields to a text file in the format:
    # field: value
    with open(outfilepath, "w") as f:
        f.write(f"num_proteins: {num_proteins}\n")
        f.write(f"mean_protein_length: {mean_protein_length}\n")
        f.write(f"median_protein_length: {median_protein_length}\n")
        f.write(f"longest_protein: {longest_protein}\n")
        f.write(f"smallest_protein: {smallest_protein}\n")
        f.write(f"from_rbh: {from_rbh}\n")

def gen_genome_stats(genomefilepath, outfilepath):
    """
    This generates statistics about a genome assembly.
    The stats that are output are:
     - the number of scaffolds
     - the GC content
     - the genome size
     - the median scaffold length
     - the mean scaffold length
     - scaffold N50
     - longest scaffold
     - smallest scaffold
     - percent Ns

    The parameters are:
      - genomefilepath: The path to the genome fasta file
      - outfilepath: The path to the output file
    The output:
      - a key: value text file with the above fields.
    """
    entries= []
    for record in fasta.parse(genomefilepath):
        entries.append({"scafname": record.id,
                        "scaflen" : len(record.seq),
                        "gc" : (record.seq.count("G") + record.seq.count("C")) / len(record.seq),
                        "Ns" : record.seq.count("N"),
                        # gaps are the number of sequential Ns of length 10 or more
                        "num_gaps": len(record.seq.upper().split("NNNNNNNNNN")) - 1})

    # make a dataframe from the entries
    df = pd.DataFrame(entries)

    num_scaffolds = len(df)
    GC_content = df["gc"].mean()
    genome_size = df["scaflen"].sum()
    median_scaffold_length = df["scaflen"].median()
    mean_scaffold_length = df["scaflen"].mean()
    scaffold_N50 = df["scaflen"].sort_values(ascending=False).cumsum().searchsorted(genome_size/2)
    longest_scaffold  = df["scaflen"].max()
    smallest_scaffold = df["scaflen"].min()
    fraction_Ns = df["Ns"].sum() / genome_size
    number_of_gaps = df["num_gaps"].sum()
    # print all of these fields to a text file in the format:
    # field: value
    with open(outfilepath, "w") as f:
        f.write(f"num_scaffolds: {num_scaffolds}\n")
        f.write(f"GC_content: {GC_content}\n")
        f.write(f"genome_size: {genome_size}\n")
        f.write(f"median_scaffold_length: {median_scaffold_length}\n")
        f.write(f"mean_scaffold_length: {mean_scaffold_length}\n")
        f.write(f"scaffold_N50: {scaffold_N50}\n")
        f.write(f"longest_scaffold: {longest_scaffold}\n")
        f.write(f"smallest_scaffold: {smallest_scaffold}\n")
        f.write(f"fraction_Ns: {fraction_Ns}\n")
        f.write(f"number_of_gaps: {number_of_gaps}\n")

def stats_filepath_to_dict(stats_filepath):
    """
    This reads in a stats file and returns a dictionary of the stats
    """
    entries = {}
    with open(stats_filepath, "r") as f:
        for line in f:
            line = line.strip()
            if line:
                key, value = line.split(": ")
                # check if the value can be cast to a float
                if value.replace(".","").isdigit():
                    if "." in value:
                        value = float(value)
                    else:
                        value = int(value)
                else:
                    # check if it is a boolean
                    if value in ["True", "False"]:
                        # we have to do this because bool("False") evaluates to True
                        if value == "True":
                            value = True
                        elif value == "False":
                            value = False
                entries[key] = value
    return entries
