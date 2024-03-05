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
        # We need to calculate
        #  - the number of scaffolds
        #  - the GC content
        #  - the genome size
        #  - the median scaffold length
        #  - the mean scaffold length
        #  - scaffold N50
        #  - longest scaffold
        #  - smallest scaffold
        #  - percent Ns
        entries= []
        for record in fasta.parse(input.genome):
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
        with open(output.results, "w") as f:
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

rule annotation_stats:
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
        # We need to calculate:
        #  - the number of proteins
        #  - the mean protein length
        #  - the median protein length
        #  - the longest protein
        #  - the smallest protein
        #  - whether the proteins are from a real annotation or from the RBH entries
        # read in the ALG_rbh file as a pandas df
        df = pd.read_csv(input.alg_rbh, sep="\t")
        rbh_names = list(df["rbh"])
        # read in the proteins. Make a list of putative rbh proteins. Get the other stats.
        entries = []
        for record in fasta.parse(input.proteins):
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
        with open(output.results, "w") as f:
            f.write(f"num_proteins: {num_proteins}\n")
            f.write(f"mean_protein_length: {mean_protein_length}\n")
            f.write(f"median_protein_length: {median_protein_length}\n")
            f.write(f"longest_protein: {longest_protein}\n")
            f.write(f"smallest_protein: {smallest_protein}\n")
            f.write(f"from_rbh: {from_rbh}\n")

rule rbhstats:
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
        # Before we calculate anything, we need to calculate:
        #  - the number of proteins in the rbh file
        #  - the number of gene groups in the rbh file
        #  - the number of genes for each gene group
        # For the final statistics, we need to calculate:
        #  - frac_ologs:           The fraction of genes of ANY ALG that are present at all in the rbh file. len(rbhdf) / total_genes_ALGs
        #  - frac_ologs_sig:       The fraction of genes of ANY ALG that are significantly on any chromosome, as defined by whole_FET
        #  - frac_ologs_single:    The fraction of genes of ANY ALG that are significantly on the largest chromosome, as defined by whole_FET
        #  - frac_ologs_{ALGNAME}: The fraction of genes of INDIVIDUAL ALGs that are significantly on any chromosome
        ALGrbhdf = pd.read_csv(input.alg_rbh, sep="\t")
        total_genes_ALGs = len(ALGrbhdf)
        genes_per_ALG    = ALGrbhdf.groupby("gene_group").size().to_dict()

        # now parse the rbh file
        rbhdf = pd.read_csv(input.rbh, sep="\t")
        sigdf = rbhdf[rbhdf["whole_FET"] < 0.05]
        # make sure that config["ALG_name"] has the appropriate columns
        for ending in ["_scaf", "_pos", "_gene"]:
            if not f"{config['ALG_name']}{ending}" in rbhdf.columns:
                raise ValueError(f"Column {config['ALG_name']}{ending} not in rbhdf.columns")
        # get the sample scaf column
        sample_loc_col = [col for col in rbhdf.columns if (col.endswith("_scaf")) and (config["ALG_name"] not in col)][0]
        frac_ologs = len(rbhdf)/total_genes_ALGs
        frac_ologs_sig = len(rbhdf[rbhdf["whole_FET"] < 0.05])/total_genes_ALGs
        # groupby the gene_group, then get the rows with the most frequent sample_loc_col value
        entries = []
        for gene_group, groupdf in sigdf.groupby("gene_group"):
            max_single = 0 if len(groupdf) == 0 else groupdf[sample_loc_col].value_counts().max()
            entries.append({"gene_group": gene_group,
                            "max_single": max_single,
                            "genes_in_group": genes_per_ALG[gene_group]})
        # frac_ologs_single is the sum o
        if len(entries) == 0:
            frac_ologs_single = float(0)
        else:
            frac_ologs_single = pd.DataFrame(entries)["max_single"].sum() / total_genes_ALGs
        # print all of these fields to a text file in the format:
        # field: value
        with open(output.results, "w") as f:
            f.write(f"frac_ologs: {frac_ologs}\n")
            f.write(f"frac_ologs_sig: {frac_ologs_sig}\n")
            f.write(f"frac_ologs_single: {frac_ologs_single}\n")
            for ALG in genes_per_ALG:
                f.write(f"frac_ologs_{ALG}: {len(sigdf[sigdf['gene_group'] == ALG])/genes_per_ALG[ALG]}\n")

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
            mydict.update(stats_filepath_to_dict(prot_stats_file))
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
            mydict.update(stats_filepath_to_dict(rbh_stats_file))
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
            mydict.update(stats_filepath_to_dict(genome_stats_file))
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