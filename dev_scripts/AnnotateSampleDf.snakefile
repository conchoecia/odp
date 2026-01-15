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

import os
import pandas as pd
import sys

# get the path of this script
snakefile_path = os.path.dirname(os.path.realpath(workflow.snakefile))

# Add paths for custom modules
# Add the parent directory to access 'source' module  
sys.path.insert(1, os.path.join(snakefile_path, ".."))
# Add dependencies path for fasta-parser
dependencies_path = os.path.join(snakefile_path, "../dependencies/fasta-parser")
sys.path.insert(1, dependencies_path)

import fasta
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
        expand(ofix + "/{dfname}.supplemented.df",
            dfname=dfname_to_filepath.keys()),
        expand(ofix + "/{dfname}.supplemented.features.pdf",
            dfname=dfname_to_filepath.keys()),
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


def genstats_get_mem_mb(wildcards, attempt):
    """
    The amount of RAM needed depends on the size of the genome file.
    """
    import os
    
    # Get the genome file path
    genome_file = config["species"].get(wildcards.sample, {}).get("genome")
    
    # Base memory allocation per attempt (in MB)
    base_mem = {
        1: 500,    # Start with 500MB
        2: 1000,   # 1GB
        3: 2500,   # 2.5GB
        4: 5000,   # 5GB
        5: 10000,  # 10GB
        6: 20000,  # 20GB
    }
    
    mem_mb = base_mem.get(attempt, 20000)
    
    # If we can check the file size, scale accordingly
    if genome_file and os.path.exists(genome_file):
        file_size_mb = os.path.getsize(genome_file) / (1024 * 1024)
        # Rule of thumb: allocate ~20x the file size in RAM for fasta parsing
        estimated_mem = int(file_size_mb * 20)
        # Use the larger of our estimate or base memory
        mem_mb = max(mem_mb, estimated_mem)
    
    return mem_mb

def genstats_get_runtime(wildcards, attempt):
    """
    The amount of runtime needed depends on the genome file size.
    """
    import os
    
    genome_file = config["species"].get(wildcards.sample, {}).get("genome")
    
    base_runtime = {
        1: 5,
        2: 10,
        3: 20,
        4: 40,
        5: 80,
        6: 160,
    }
    
    runtime = base_runtime.get(attempt, 160)
    
    if genome_file and os.path.exists(genome_file):
        file_size_mb = os.path.getsize(genome_file) / (1024 * 1024)
        # Rule of thumb: 1 minute per 50MB of genome
        estimated_runtime = int(file_size_mb / 50)
        runtime = max(5, max(runtime, estimated_runtime))
    
    return runtime

def annotation_stats_get_mem_mb(wildcards, attempt):
    """
    The amount of RAM needed depends on the size of the protein file.
    """
    import os
    
    # Get the protein file path
    protein_file = config["species"].get(wildcards.sample, {}).get("proteins")
    
    # Base memory allocation per attempt (in MB)
    base_mem = {
        1: 1000,   # Start with 1GB
        2: 2500,   # 2.5GB
        3: 5000,   # 5GB
        4: 10000,  # 10GB
        5: 20000,  # 20GB
        6: 40000,  # 40GB
    }
    
    mem_mb = base_mem.get(attempt, 40000)
    
    # If we can check the file size, scale accordingly
    if protein_file and os.path.exists(protein_file):
        file_size_mb = os.path.getsize(protein_file) / (1024 * 1024)
        # Rule of thumb: allocate ~30x the file size in RAM
        estimated_mem = int(file_size_mb * 30)
        # Use the larger of our estimate or base memory
        mem_mb = max(mem_mb, estimated_mem)
    
    return mem_mb

def annotation_stats_get_runtime(wildcards, attempt):
    """
    The amount of runtime needed depends on the protein file size.
    """
    import os
    
    protein_file = config["species"].get(wildcards.sample, {}).get("proteins")
    
    base_runtime = {
        1: 5,
        2: 10,
        3: 20,
        4: 40,
        5: 80,
        6: 160,
    }
    
    runtime = base_runtime.get(attempt, 160)
    
    if protein_file and os.path.exists(protein_file):
        file_size_mb = os.path.getsize(protein_file) / (1024 * 1024)
        # Rule of thumb: 1 minute per 20MB of protein file
        estimated_runtime = int(file_size_mb / 20)
        runtime = max(5, max(runtime, estimated_runtime))
    
    return runtime

def rbhstats_get_mem_mb(wildcards, attempt):
    """
    The amount of RAM needed depends on the size of the RBH file.
    RBH files can be quite large, so we allocate more RAM by default.
    We also check the file size to scale appropriately.
    """
    import os
    
    # Get the RBH file path
    rbh_file = config["species_to_rbh"].get(wildcards.sample)
    
    # Base memory allocation per attempt (in MB)
    base_mem = {
        1: 2000,   # Start with 2GB instead of 250MB
        2: 5000,   # 5GB
        3: 10000,  # 10GB
        4: 20000,  # 20GB
        5: 40000,  # 40GB
        6: 80000,  # 80GB
    }
    
    mem_mb = base_mem.get(attempt, 80000)
    
    # If we can check the file size, scale accordingly
    if rbh_file and os.path.exists(rbh_file):
        file_size_mb = os.path.getsize(rbh_file) / (1024 * 1024)
        # Rule of thumb: allocate ~50x the file size in RAM
        estimated_mem = int(file_size_mb * 50)
        # Use the larger of our estimate or base memory
        mem_mb = max(mem_mb, estimated_mem)
    
    return mem_mb

def rbhstats_get_runtime(wildcards, attempt):
    """
    The amount of runtime needed depends on the file size and attempt number.
    """
    import os
    
    # Get the RBH file path
    rbh_file = config["species_to_rbh"].get(wildcards.sample)
    
    # Base runtime allocation per attempt (in minutes)
    base_runtime = {
        1: 10,   # Start with 10 minutes instead of 2
        2: 20,
        3: 40,
        4: 80,
        5: 160,
        6: 320,
    }
    
    runtime = base_runtime.get(attempt, 320)
    
    # If we can check the file size, scale accordingly
    if rbh_file and os.path.exists(rbh_file):
        file_size_mb = os.path.getsize(rbh_file) / (1024 * 1024)
        # Rule of thumb: 1 minute per 10MB of file
        estimated_runtime = int(file_size_mb / 10)
        # Use the larger of our estimate or base runtime, with a minimum of 5 minutes
        runtime = max(5, max(runtime, estimated_runtime))
    
    return runtime

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
    log:
        ofix + "/logs/genome_stats/{sample}.genome_stats.log"
    benchmark:
        ofix + "/benchmarks/genome_stats/{sample}.genome_stats.tsv"
    threads: 1
    retries: 6
    resources:
        mem_mb = genstats_get_mem_mb,
        runtime = genstats_get_runtime
    run:
        import os
        import sys
        from datetime import datetime
        
        os.makedirs(os.path.dirname(log[0]), exist_ok=True)
        
        genome_size_mb = os.path.getsize(input.genome) / (1024 * 1024)
        
        with open(log[0], 'w') as logfile:
            logfile.write(f"Starting genstats for sample: {wildcards.sample}\n")
            logfile.write(f"Timestamp: {datetime.now().isoformat()}\n")
            logfile.write(f"Genome file: {input.genome}\n")
            logfile.write(f"Genome file size: {genome_size_mb:.2f} MB\n")
            logfile.write(f"Allocated memory: {resources.mem_mb} MB\n")
            logfile.write(f"Allocated runtime: {resources.runtime} minutes\n")
            logfile.write("-" * 80 + "\n")
            
            old_stdout = sys.stdout
            old_stderr = sys.stderr
            sys.stdout = logfile
            sys.stderr = logfile
            
            try:
                asd.gen_genome_stats(input.genome, output.results)
                logfile.write("\n" + "-" * 80 + "\n")
                logfile.write(f"Completed successfully at {datetime.now().isoformat()}\n")
            except Exception as e:
                logfile.write("\n" + "-" * 80 + "\n")
                logfile.write(f"ERROR: {str(e)}\n")
                logfile.write(f"Failed at {datetime.now().isoformat()}\n")
                raise
            finally:
                sys.stdout = old_stdout
                sys.stderr = old_stderr

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
    log:
        ofix + "/logs/protein_stats/{sample}.protein_stats.log"
    benchmark:
        ofix + "/benchmarks/protein_stats/{sample}.protein_stats.tsv"
    threads: 1
    retries: 6
    resources:
        mem_mb = annotation_stats_get_mem_mb,
        runtime = annotation_stats_get_runtime
    run:
        import os
        import sys
        from datetime import datetime
        
        os.makedirs(os.path.dirname(log[0]), exist_ok=True)
        
        protein_size_mb = os.path.getsize(input.proteins) / (1024 * 1024)
        alg_size_mb = os.path.getsize(input.alg_rbh) / (1024 * 1024)
        
        with open(log[0], 'w') as logfile:
            logfile.write(f"Starting annotation_stats for sample: {wildcards.sample}\n")
            logfile.write(f"Timestamp: {datetime.now().isoformat()}\n")
            logfile.write(f"Protein file: {input.proteins}\n")
            logfile.write(f"Protein file size: {protein_size_mb:.2f} MB\n")
            logfile.write(f"ALG RBH file: {input.alg_rbh}\n")
            logfile.write(f"ALG RBH file size: {alg_size_mb:.2f} MB\n")
            logfile.write(f"Allocated memory: {resources.mem_mb} MB\n")
            logfile.write(f"Allocated runtime: {resources.runtime} minutes\n")
            logfile.write("-" * 80 + "\n")
            
            old_stdout = sys.stdout
            old_stderr = sys.stderr
            sys.stdout = logfile
            sys.stderr = logfile
            
            try:
                asd.gen_annotation_stats(input.proteins, input.alg_rbh, output.results)
                logfile.write("\n" + "-" * 80 + "\n")
                logfile.write(f"Completed successfully at {datetime.now().isoformat()}\n")
            except Exception as e:
                logfile.write("\n" + "-" * 80 + "\n")
                logfile.write(f"ERROR: {str(e)}\n")
                logfile.write(f"Failed at {datetime.now().isoformat()}\n")
                raise
            finally:
                sys.stdout = old_stdout
                sys.stderr = old_stderr

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
    log:
        ofix + "/logs/rbh_stats/{sample}.rbh_stats.log"
    benchmark:
        ofix + "/benchmarks/rbh_stats/{sample}.rbh_stats.tsv"
    threads: 1
    retries: 6
    resources:
        mem_mb  = rbhstats_get_mem_mb,
        runtime = rbhstats_get_runtime
    run:
        import os
        import sys
        from datetime import datetime
        
        # Create log directory if it doesn't exist
        os.makedirs(os.path.dirname(log[0]), exist_ok=True)
        
        # Log file size and resource allocation
        rbh_size_mb = os.path.getsize(input.rbh) / (1024 * 1024)
        alg_size_mb = os.path.getsize(input.alg_rbh) / (1024 * 1024)
        
        with open(log[0], 'w') as logfile:
            logfile.write(f"Starting rbhstats for sample: {wildcards.sample}\n")
            logfile.write(f"Timestamp: {datetime.now().isoformat()}\n")
            logfile.write(f"RBH file: {input.rbh}\n")
            logfile.write(f"RBH file size: {rbh_size_mb:.2f} MB\n")
            logfile.write(f"ALG RBH file: {input.alg_rbh}\n")
            logfile.write(f"ALG RBH file size: {alg_size_mb:.2f} MB\n")
            logfile.write(f"Allocated memory: {resources.mem_mb} MB\n")
            logfile.write(f"Allocated runtime: {resources.runtime} minutes\n")
            logfile.write(f"Threads: {threads}\n")
            logfile.write("-" * 80 + "\n")
            
            # Redirect stdout and stderr to log file
            old_stdout = sys.stdout
            old_stderr = sys.stderr
            sys.stdout = logfile
            sys.stderr = logfile
            
            try:
                asd.gen_rbh_stats(input.rbh, input.alg_rbh, config["ALG_name"], output.results)
                logfile.write("\n" + "-" * 80 + "\n")
                logfile.write(f"Completed successfully at {datetime.now().isoformat()}\n")
            except Exception as e:
                logfile.write("\n" + "-" * 80 + "\n")
                logfile.write(f"ERROR: {str(e)}\n")
                logfile.write(f"Failed at {datetime.now().isoformat()}\n")
                raise
            finally:
                sys.stdout = old_stdout
                sys.stderr = old_stderr

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
        mem_mb  = 8000,
        runtime = 10
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
        pdf = ofix + "/{dfname}.supplemented.features.pdf"
    threads: 1
    resources:
        mem_mb  = 4000,
        runtime = 10
    params:
        prefix = lambda wildcards, output: output.pdf.replace('.features.pdf', '')
    shell:
        """
        python {input.plotdfs} --plot_features -f {input.df} -p {params.prefix} --threecolor --genome-min-bp 100000000 --genome-max-bp 5000000000 --num-cols 7
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

rule collate_rbhstats_logs:
    """
    Collates the rbhstats log files into a single TSV for analysis.
    Extracts: sample, RBH file size, ALG RBH file size, allocated memory, allocated runtime
    """
    input:
        logs = expand(ofix + "/logs/rbh_stats/{sample}.rbh_stats.log", sample=config["species"].keys()),
        benchmarks = expand(ofix + "/benchmarks/rbh_stats/{sample}.rbh_stats.tsv", sample=config["species"].keys())
    output:
        tsv = ofix + "/logs/rbhstats_resource_usage.tsv"
    threads: 1
    resources:
        mem_mb = 2000,
        runtime = 30
    run:
        import re
        import os
        import pandas as pd
        
        results = []
        
        for log_file in input.logs:
            sample = os.path.basename(log_file).replace('.rbh_stats.log', '')
            
            # Parse log file
            log_data = {
                'sample': sample,
                'rbh_size_mb': None,
                'alg_rbh_size_mb': None,
                'allocated_mem_mb': None,
                'allocated_runtime_min': None,
                'actual_runtime_sec': None,
                'actual_mem_mb': None,
                'success': False
            }
            
            try:
                with open(log_file, 'r') as f:
                    content = f.read()
                    
                    # Extract file sizes
                    rbh_match = re.search(r'RBH file size: ([\d.]+) MB', content)
                    if rbh_match:
                        log_data['rbh_size_mb'] = float(rbh_match.group(1))
                    
                    alg_match = re.search(r'ALG RBH file size: ([\d.]+) MB', content)
                    if alg_match:
                        log_data['alg_rbh_size_mb'] = float(alg_match.group(1))
                    
                    # Extract allocated resources
                    mem_match = re.search(r'Allocated memory: (\d+) MB', content)
                    if mem_match:
                        log_data['allocated_mem_mb'] = int(mem_match.group(1))
                    
                    runtime_match = re.search(r'Allocated runtime: (\d+) minutes', content)
                    if runtime_match:
                        log_data['allocated_runtime_min'] = int(runtime_match.group(1))
                    
                    # Check success
                    if 'Completed successfully' in content:
                        log_data['success'] = True
            except Exception as e:
                print(f"Warning: Could not parse log file {log_file}: {e}")
            
            # Parse benchmark file for actual usage
            benchmark_file = log_file.replace('/logs/rbh_stats/', '/benchmarks/rbh_stats/').replace('.log', '.tsv')
            try:
                if os.path.exists(benchmark_file):
                    df_bench = pd.read_csv(benchmark_file, sep='\t')
                    if not df_bench.empty:
                        log_data['actual_runtime_sec'] = df_bench['s'].iloc[0]
                        log_data['actual_mem_mb'] = df_bench['max_rss'].iloc[0] / 1024  # Convert KB to MB
            except Exception as e:
                print(f"Warning: Could not parse benchmark file {benchmark_file}: {e}")
            
            results.append(log_data)
        
        # Create DataFrame and save
        df = pd.DataFrame(results)
        df = df.sort_values('sample')
        df.to_csv(output.tsv, sep='\t', index=False)
        
        # Print summary statistics
        print(f"\nRBHStats Resource Usage Summary:")
        print(f"Total samples: {len(df)}")
        print(f"Successful: {df['success'].sum()}")
        print(f"Failed: {(~df['success']).sum()}")
        if df['actual_mem_mb'].notna().any():
            print(f"\nMemory Usage (MB):")
            print(f"  Mean: {df['actual_mem_mb'].mean():.2f}")
            print(f"  Median: {df['actual_mem_mb'].median():.2f}")
            print(f"  Max: {df['actual_mem_mb'].max():.2f}")
        if df['actual_runtime_sec'].notna().any():
            print(f"\nRuntime (seconds):")
            print(f"  Mean: {df['actual_runtime_sec'].mean():.2f}")
            print(f"  Median: {df['actual_runtime_sec'].median():.2f}")
            print(f"  Max: {df['actual_runtime_sec'].max():.2f}")

rule collate_genstats_logs:
    """
    Collates the genstats log files into a single TSV for analysis.
    Extracts: sample, genome file size, allocated memory, allocated runtime
    """
    input:
        logs = expand(ofix + "/logs/genome_stats/{sample}.genome_stats.log", sample=config["species"].keys()),
        benchmarks = expand(ofix + "/benchmarks/genome_stats/{sample}.genome_stats.tsv", sample=config["species"].keys())
    output:
        tsv = ofix + "/logs/genstats_resource_usage.tsv"
    threads: 1
    resources:
        mem_mb = 2000,
        runtime = 30
    run:
        import re
        import os
        import pandas as pd
        
        results = []
        
        for log_file in input.logs:
            sample = os.path.basename(log_file).replace('.genome_stats.log', '')
            
            log_data = {
                'sample': sample,
                'genome_size_mb': None,
                'allocated_mem_mb': None,
                'allocated_runtime_min': None,
                'actual_runtime_sec': None,
                'actual_mem_mb': None,
                'success': False
            }
            
            try:
                with open(log_file, 'r') as f:
                    content = f.read()
                    
                    # Extract file size
                    size_match = re.search(r'Genome file size: ([\d.]+) MB', content)
                    if size_match:
                        log_data['genome_size_mb'] = float(size_match.group(1))
                    
                    # Extract allocated resources
                    mem_match = re.search(r'Allocated memory: (\d+) MB', content)
                    if mem_match:
                        log_data['allocated_mem_mb'] = int(mem_match.group(1))
                    
                    runtime_match = re.search(r'Allocated runtime: (\d+) minutes', content)
                    if runtime_match:
                        log_data['allocated_runtime_min'] = int(runtime_match.group(1))
                    
                    # Check success
                    if 'Completed successfully' in content:
                        log_data['success'] = True
            except Exception as e:
                print(f"Warning: Could not parse log file {log_file}: {e}")
            
            # Parse benchmark file
            benchmark_file = log_file.replace('/logs/genome_stats/', '/benchmarks/genome_stats/').replace('.log', '.tsv')
            try:
                if os.path.exists(benchmark_file):
                    df_bench = pd.read_csv(benchmark_file, sep='\t')
                    if not df_bench.empty:
                        log_data['actual_runtime_sec'] = df_bench['s'].iloc[0]
                        log_data['actual_mem_mb'] = df_bench['max_rss'].iloc[0] / 1024
            except Exception as e:
                print(f"Warning: Could not parse benchmark file {benchmark_file}: {e}")
            
            results.append(log_data)
        
        df = pd.DataFrame(results)
        df = df.sort_values('sample')
        df.to_csv(output.tsv, sep='\t', index=False)
        
        print(f"\nGenStats Resource Usage Summary:")
        print(f"Total samples: {len(df)}")
        print(f"Successful: {df['success'].sum()}")
        print(f"Failed: {(~df['success']).sum()}")
        if df['actual_mem_mb'].notna().any():
            print(f"\nMemory Usage (MB):")
            print(f"  Mean: {df['actual_mem_mb'].mean():.2f}")
            print(f"  Median: {df['actual_mem_mb'].median():.2f}")
            print(f"  Max: {df['actual_mem_mb'].max():.2f}")
        if df['actual_runtime_sec'].notna().any():
            print(f"\nRuntime (seconds):")
            print(f"  Mean: {df['actual_runtime_sec'].mean():.2f}")
            print(f"  Median: {df['actual_runtime_sec'].median():.2f}")
            print(f"  Max: {df['actual_runtime_sec'].max():.2f}")

rule collate_annotation_stats_logs:
    """
    Collates the annotation_stats log files into a single TSV for analysis.
    Extracts: sample, protein file size, ALG RBH file size, allocated memory, allocated runtime
    """
    input:
        logs = expand(ofix + "/logs/protein_stats/{sample}.protein_stats.log", sample=config["species"].keys()),
        benchmarks = expand(ofix + "/benchmarks/protein_stats/{sample}.protein_stats.tsv", sample=config["species"].keys())
    output:
        tsv = ofix + "/logs/annotation_stats_resource_usage.tsv"
    threads: 1
    resources:
        mem_mb = 2000,
        runtime = 30
    run:
        import re
        import os
        import pandas as pd
        
        results = []
        
        for log_file in input.logs:
            sample = os.path.basename(log_file).replace('.protein_stats.log', '')
            
            log_data = {
                'sample': sample,
                'protein_size_mb': None,
                'alg_rbh_size_mb': None,
                'allocated_mem_mb': None,
                'allocated_runtime_min': None,
                'actual_runtime_sec': None,
                'actual_mem_mb': None,
                'success': False
            }
            
            try:
                with open(log_file, 'r') as f:
                    content = f.read()
                    
                    # Extract file sizes
                    protein_match = re.search(r'Protein file size: ([\d.]+) MB', content)
                    if protein_match:
                        log_data['protein_size_mb'] = float(protein_match.group(1))
                    
                    alg_match = re.search(r'ALG RBH file size: ([\d.]+) MB', content)
                    if alg_match:
                        log_data['alg_rbh_size_mb'] = float(alg_match.group(1))
                    
                    # Extract allocated resources
                    mem_match = re.search(r'Allocated memory: (\d+) MB', content)
                    if mem_match:
                        log_data['allocated_mem_mb'] = int(mem_match.group(1))
                    
                    runtime_match = re.search(r'Allocated runtime: (\d+) minutes', content)
                    if runtime_match:
                        log_data['allocated_runtime_min'] = int(runtime_match.group(1))
                    
                    # Check success
                    if 'Completed successfully' in content:
                        log_data['success'] = True
            except Exception as e:
                print(f"Warning: Could not parse log file {log_file}: {e}")
            
            # Parse benchmark file
            benchmark_file = log_file.replace('/logs/protein_stats/', '/benchmarks/protein_stats/').replace('.log', '.tsv')
            try:
                if os.path.exists(benchmark_file):
                    df_bench = pd.read_csv(benchmark_file, sep='\t')
                    if not df_bench.empty:
                        log_data['actual_runtime_sec'] = df_bench['s'].iloc[0]
                        log_data['actual_mem_mb'] = df_bench['max_rss'].iloc[0] / 1024
            except Exception as e:
                print(f"Warning: Could not parse benchmark file {benchmark_file}: {e}")
            
            results.append(log_data)
        
        df = pd.DataFrame(results)
        df = df.sort_values('sample')
        df.to_csv(output.tsv, sep='\t', index=False)
        
        print(f"\nAnnotation Stats Resource Usage Summary:")
        print(f"Total samples: {len(df)}")
        print(f"Successful: {df['success'].sum()}")
        print(f"Failed: {(~df['success']).sum()}")
        if df['actual_mem_mb'].notna().any():
            print(f"\nMemory Usage (MB):")
            print(f"  Mean: {df['actual_mem_mb'].mean():.2f}")
            print(f"  Median: {df['actual_mem_mb'].median():.2f}")
            print(f"  Max: {df['actual_mem_mb'].max():.2f}")
        if df['actual_runtime_sec'].notna().any():
            print(f"\nRuntime (seconds):")
            print(f"  Mean: {df['actual_runtime_sec'].mean():.2f}")
            print(f"  Median: {df['actual_runtime_sec'].median():.2f}")
            print(f"  Max: {df['actual_runtime_sec'].max():.2f}")
