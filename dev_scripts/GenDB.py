#!/usr/bin/env python

"""
This directory 
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