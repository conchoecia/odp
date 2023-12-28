#!/usr/bin/env python

"""
This python program contains functions that are used to download genomes from NCBI.
  - These tools are common to the scripts:
    - GenDB_build_db_unannotated_chr.snakefile
    - GenDB_build_db_unannotated_nonchr.snakefile
    - GenDB_build_db_annotated_chr.snakefile
    - GenDB_build_db_annotated_nonchr.snakefile
"""

import os
import pandas as pd

def opening_logic_GenDB_build_db(config, chr_scale = None, annotated = None):
    """
    This is common logic for all of the GenDB_build_db_*.snakefile scripts.
    Basically just pertains to parsing which files will go into this category.

    Takes the config object as input, and returns it modified.
    The user must provide boolean values (True, False) for chr_scale and annotated.
    """
    # Do some logic to see if the user has procided enough information for us to analyse the genomes
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
        latest_accesions = return_latest_accession_tsvs(config["directory"])
        # This is the output of the above function:
        #return_dict = {"annotated_chr":      os.path.join(directory_path, most_recent_annotated_chr_file),
        #               "annotated_nonchr":   os.path.join(directory_path, most_recent_annotated_nonchr_file),
        #               "unannotated_chr":    os.path.join(directory_path, most_recent_unannotated_chr_file),
        #               "unannotated_nonchr": os.path.join(directory_path, most_recent_unannotated_nonchr_file)
        #               }
        config["annotated_genome_chr_tsv"]       = latest_accesions["annotated_chr"]
        config["annotated_genome_nonchr_tsv"]    = latest_accesions["annotated_nonchr"]
        config["unannotated_genome_chr_tsv"]     = latest_accesions["unannotated_chr"]
        config["unannotated_genome_nonchr_tsv"]  = latest_accesions["unannotated_nonchr"]

        # Which file we use depends on whether we are doing annotated or unannotated genomes,
        #  and whether we are doing chromosome or non-chromosome scale genomes.
        chosen_file = None
        if annotated:
            if chr_scale:
                chosen_file = config["annotated_genome_chr_tsv"]
            else:
                chosen_file = config["annotated_genome_nonchr_tsv"]
        else:
            if chr_scale:
                chosen_file = config["unannotated_genome_chr_tsv"]
            else:
                chosen_file = config["unannotated_genome_nonchr_tsv"]
        # Use the chosen file, and add the entries to the config file so we can download them or not
        config["assemAnn"] = determine_genome_accessions(chosen_file)

        # TODO This line should probably get its own function later
        ## get the list of GCAs to ignore in case we need to remove any
        #ignore_list_path = os.path.join(snakefile_path, "assembly_ignore_list.txt")
        #with open(ignore_list_path, "r") as f:
        #    for line in f:
        #        line = line.strip()
        #        if line:
        #            if line in config["assemAnn"]:
        #                config["assemAnn"].remove(line)

    elif "accession_tsvs" in config:
        # ensure that the user also hasn't specified the directory
        if "directory" in config:
            raise IOerror("You cannot provide both a directory of tsv files ('directory') and specify the exact path of the TSVs ('accession_tsvs'). Read the config file.")
        # we haven't implemented this yet. I haven't found a use case where I would want to specifially pick a file path rather than just get the most recent one.
        raise NotImplementedError("We haven't implemented this yet. I haven't found a use case where I would want to specifially pick a file path rather than just get the most recent one.")

    return config


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

    Returns the full path to the most recent annotated and unannotated accession TSV files.
    """
    # check that the directory exists
    if not os.path.isdir(directory_path):
        raise IOerror("The directory of TSV files you provided does not exist. {}".format(directory_path))
    # get a list of files in this directory that start with ["annotated", "unannotated"]
    files = os.listdir(directory_path)
    annotated_chr_files      = [f for f in files if f.startswith("annotated_genomes_chr")      and f.endswith(".tsv") and contains_date(f)]
    annotated_nonchr_files   = [f for f in files if f.startswith("annotated_genomes_nonchr")   and f.endswith(".tsv") and contains_date(f)]
    unannotated_chr_files    = [f for f in files if f.startswith("unannotated_genomes_chr")    and f.endswith(".tsv") and contains_date(f)]
    unannotated_nonchr_files = [f for f in files if f.startswith("unannotated_genomes_nonchr") and f.endswith(".tsv") and contains_date(f)]
    # sort the files by date, using YYYYMMDDHHMM as the sorting key
    annotated_chr_files.sort(        key=lambda x: int(x.split("_")[-1].split(".")[0]), reverse=True)
    annotated_nonchr_files.sort(     key=lambda x: int(x.split("_")[-1].split(".")[0]), reverse=True)
    unannotated_chr_files.sort(      key=lambda x: int(x.split("_")[-1].split(".")[0]), reverse=True)
    unannotated_nonchr_files.sort(   key=lambda x: int(x.split("_")[-1].split(".")[0]), reverse=True)
    # get the most recent annotated file
    most_recent_annotated_chr_file      = annotated_chr_files[0] if len(annotated_chr_files) > 0 else "None.txt"
    most_recent_annotated_nonchr_file   = annotated_nonchr_files[0] if len(annotated_nonchr_files) > 0 else "None.txt"
    most_recent_unannotated_chr_file    = unannotated_chr_files[0]  if len(unannotated_chr_files) > 0 else "None.txt"
    most_recent_unannotated_nonchr_file = unannotated_nonchr_files[0]  if len(unannotated_nonchr_files) > 0 else "None.txt"
    return_dict = {"annotated_chr":      os.path.join(directory_path, most_recent_annotated_chr_file),
                   "annotated_nonchr":   os.path.join(directory_path, most_recent_annotated_nonchr_file),
                   "unannotated_chr":    os.path.join(directory_path, most_recent_unannotated_chr_file),
                   "unannotated_nonchr": os.path.join(directory_path, most_recent_unannotated_nonchr_file)
                   }
    return return_dict

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
