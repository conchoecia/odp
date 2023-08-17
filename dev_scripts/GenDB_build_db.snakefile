"""
This takes the list of annotated and unannotated genomes and prepares a database from them for ODP.

The program requires that either a directory of the annotated and unannotated genome lists be provided.

Otherwise the user has to specify specific paths to those tsv files.
"""

configfile: "config.yaml"

config["tool"] = "odp_ncbi_genome_db"

wildcard_constraints:
    taxid="[0-9]+",

import os
import pandas as pd
import sys
from datetime import datetime

#if "datetime" not in config:
#    config["datetime"] = datetime.now().strftime('%Y%m%d%H%M')

# figure out where bin is because we need to use some outside tools
snakefile_path = os.path.dirname(os.path.realpath(workflow.snakefile))
bin_path = os.path.join(snakefile_path, "../bin")

def contains_date(string_to_check):
    """
    From the string in question, just see if it contains a date in the format YYYYMMDDHHMM
    """
    # split the string on the underscore
    split_string = string_to_check.replace(".tsv","").split("_")
    # check if any of the elements are a date
    string_contains_date = False
    for element in split_string:
        if element.isdigit() and len(element) == 12:
            string_contains_date = True
    return string_contains_date

def return_latest_accession_tsvs(directory_path):
    """
    Given a directory, we have to do some parsing of the file names to figure out what is the most recent accession tsv file.
    """
    # check that the directory exists
    if not os.path.isdir(directory_path):
        raise IOerror("The directory of TSV files you provided does not exist. {}".format(directory_path))
    # get a list of files in this directory that start with ["annotated", "unannotated"]
    files = os.listdir(directory_path)
    annotated_files   = [f for f in files if f.startswith("annotated") and f.endswith(".tsv") and contains_date(f)]
    unannotated_files = [f for f in files if f.startswith("unannotated") and f.endswith(".tsv") and contains_date(f)]
    # sort the files by date, using YYYYMMDDHHMM as the sorting key
    annotated_files.sort(key=lambda x: int(x.split("_")[-1].split(".")[0]),    reverse=True)
    unannotated_files.sort(key=lambda x: int(x.split("_")[-1].split(".")[0]) , reverse=True)
    # get the most recent annotated file
    most_recent_annotated_file = annotated_files[0]
    most_recent_unannotated_file = unannotated_files[0]
    return os.path.join(directory_path, most_recent_annotated_file), os.path.join(directory_path, most_recent_unannotated_file)

def determine_genome_accessions(tsv_filepath):
    """
    Reads in the TSV file that has the information about the genomes,
    and returns a data structure (TBD) of which sequences to download.

    This generic function can be used to get information both about the unannotated and annotated genomes.

    Just returns the "Assembly Accession" column of the TSV file as a list.
    This assembly accession number is all that is needed to download the genome and the annotation.
    """
    df = pd.read_csv(tsv_filepath, sep="\t")
    # strip leading and trailing whitespace from the column names because pandas can screw up sometimes 
    df.columns = df.columns.str.strip()
    return df["Assembly Accession"].tolist()

# Do some logic to see if the user has procided enough informatoin for us to analyse the genomes
if ("directory" not in config) and ("accession_tsvs" not in config): 
    raise IOerror("You must provide either a directory of the annotated and unannotated genome lists, or a list of the paths to those tsv files. Read the config file.")
if "directory" in config:
    # ensure that the user also hasn't specified the accession tsvs
    if "accession_tsvs" in config:
        raise IOerror("You cannot provide both a directory of tsv files ('directory') and specify the exact path of the TSVs ('accession_tsvs'). Read the config file.")
    # ensure that the directory exists
    if not os.path.isdir(config["directory"]):
        raise IOerror("The directory of TSV files you provided does not exist. {}".format(config["directory"]))
    # get the paths to the tsv files
    config["annotated_genome_tsv"], config["unannotated_genome_tsv"] = return_latest_accession_tsvs(config["directory"])
    # now add the entries to the config file so we can download them or not
    config["assemAnn"] = determine_genome_accessions(config["annotated_genome_tsv"])

elif "accession_tsvs" in config:
    # ensure that the user also hasn't specified the directory
    if "directory" in config:
        raise IOerror("You cannot provide both a directory of tsv files ('directory') and specify the exact path of the TSVs ('accession_tsvs'). Read the config file.")
    # we haven't implemented this yet. I haven't found a use case where I would want to specifially pick a file path rather than just get the most recent one.
    raise NotImplementedError("We haven't implemented this yet. I haven't found a use case where I would want to specifially pick a file path rather than just get the most recent one.")

# first we must load in all of the files. Only do it once
rule all:
    input:
        #expand(config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.yaml.part", assemAnn=config["assemAnn"]),
        "NCBI_odp_db.yaml"

rule download_annotated_genomes:
    """
    We have selected the annotated genomes to download. These are the easiest to handle since we don't have to annotate them ourselves.
    """
    input:
        datasets = os.path.join(bin_path, "datasets")
    output:
        assembly = temp(config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/ncbi_dataset.zip")
    params:
        outdir   = config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/"
    threads: 1
    resources:
        mem_mb = 1000
    shell:
        """
        cd {params.outdir}
        {input.datasets} download genome accession {wildcards.assemAnn} --include genome,protein,gff3,gtf
        """

rule unzip_annotated_genomes:
    """
    We just downloaded the genome data packet. Now unzip it, delete the zip file, and rename things
    """
    input:
        assembly = config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/ncbi_dataset.zip"
    output:
        genome  = config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.fasta",
        protein = config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.pep",
        gff     = config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.gff",
        gtf     = config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.gtf"
    threads: 1
    resources:
        mem_mb = 1000
    params:
        outdir   = config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/"
    shell:
        """
        TMPDIR=`pwd`
        cd {params.outdir}
        unzip -o ncbi_dataset.zip
        cd $TMPDIR
        find {params.outdir} -name "*.fna" -exec mv {{}} {output.genome} \;
        find {params.outdir} -name "*.faa" -exec mv {{}} {output.protein} \;
        find {params.outdir} -name "*.gff" -exec mv {{}} {output.gff} \;
        find {params.outdir} -name "*.gtf" -exec mv {{}} {output.gtf} \;
        """

rule prep_chrom_file_from_NCBI:
    """
    This takes the output files from the NCBI database and puts them into the file format we need for odp, clink, et cetera
    """
    input:
        genome   = config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.fasta",
        protein  = config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.pep",
        gff      = config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.gff",
        chromgen = os.path.join(snakefile_path, "..", "scripts", "NCBIgff2chrom.py")
    output:
        chrom   = config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.chrom",
        report  = config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.report.txt"
    threads: 1
    resources:
        mem_mb  = 1000
    params:
        prefix  = config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}"
    shell:
        """
        python {input.chromgen} -f {input.genome} -p {input.protein} -g {input.gff} -o {params.prefix}
        """

rule generate_assembled_config_entry:
    """
    Print out a small piece of a yaml file specifically for ODP.
    These will be gathered and concatenated later.
    """
    input:
        annotated_genomes = config["annotated_genome_tsv"],
        report            = config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.report.txt",
        genome            = config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.fasta",
        protein           = config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.pep",
        chrom             = config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.chrom",
    output:
        yaml   = config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.yaml.part",
    threads: 1
    resources:
        mem_mb  = 1000
    run:
        # load in the dataframe of the annotated genomes
        df = pd.read_csv(input.annotated_genomes, sep="\t")
        # strip leading and trailing whitespace from the column names because pandas can screw up sometimes
        df.columns = df.columns.str.strip()

        # read in the report into pandas if we can. Strip all the excess whitespace since the columns are formatted with whitespace
        reportdf = pd.read_csv(input.report, comment = "#", delim_whitespace=True)
        # get only the scaffolds that comprise more than 0.5% of the annotated proteins
        filtdf = reportdf.loc[reportdf["percent_of_proteins"] >= 0.05]
        minscaflen = filtdf["scaflen"].min() - 1000

        row = df.loc[df["Assembly Accession"] == wildcards.assemAnn]
        taxid = row["Organism Taxonomic ID"].values[0]

        # s is the output string.
        # h is headspace. Currently 2 spaces because the odp yaml files are compoased like so:
        #species:
        #  Sp_name:
        #    taxid:
        #    genus:
        #    species:
        #    assembly_accession: 
        #    proteins:
        #    chrom:
        #    genome:
        #    minscaflen:
        #    ....

        # genus
        try:
            genus = row["Organism Name"].values[0].split(" ")[0]
        except:
            genus = row["Organism Name"].values[0]
        # species
        try:
            species = row["Organism Name"].values[0].split(" ")[1]
        except:
            species = "None" 
        
        spstring = "{}{}{}".format(genus, species, taxid)
        h = "  "
        s = ""
        s += h + "{}:\n".format(spstring)
        s += h + h + "assembly_accession: {}\n".format(wildcards.assemAnn)
        s += h + h + "taxid:              {}\n".format(taxid)
        s += h + h + "genus:              {}\n".format(genus)
        s += h + h + "species:            {}\n".format(species)
        s += h + h + "proteins:           {}\n".format(os.path.abspath(input.protein))
        s += h + h + "chrom:              {}\n".format(os.path.abspath(input.chrom))
        s += h + h + "genome:             {}\n".format(os.path.abspath(input.genome))
        s += h + h + "minscaflen:         {}\n".format(minscaflen)

        # write the string to the output
        with open(output.yaml, "w") as f:
            f.write(s)

rule collate_assembled_config_entries:
    """
    Concatenate all of the yaml files together into one big file.
    """
    input:
        yaml_parts = expand(config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.yaml.part",
                            assemAnn=config["assemAnn"])
    output:
        yaml = "NCBI_odp_db.yaml"
    threads: 1
    resources:
        mem_mb  = 1000
    shell:
        """
        cat {input.yaml_parts} > {output.yaml}
        """