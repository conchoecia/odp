#!/usr/bin/env python

"""
This python program contains functions that are used to download genomes from NCBI.
  - These tools are common to the scripts:
    - GenDB_build_db_unannotated_chr.snakefile
    - GenDB_build_db.snakefile
"""

import os
import pandas as pd

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
    most_recent_annotated_chr_file      = annotated_chr_files[0]
    most_recent_annotated_nonchr_file   = annotated_nonchr_files[0]
    most_recent_unannotated_chr_file    = unannotated_chr_files[0]
    most_recent_unannotated_nonchr_file = unannotated_nonchr_files[0]
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
