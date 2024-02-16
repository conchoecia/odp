#!/usr/bin/env python
"""
Program  : rbh_tools.py
Language : python
Date     : 2024-02-07
Author   : Darrin T. Schultz
Email    : darrin.schultz@univie.ac.at
Github   : https://github.com/conchoecia/odp
Support  : For issues or questions, please search if the topic has been discussed already
           on github and open a new issue if not: https://github.com/conchoecia/odp/issues
License  : GNU GENERAL PUBLIC LICENSE, Version 3, 29 June 2007. See the LICENSE file.
Citation : If you use this software for your scientific publication, please cite:
           Schultz, DT; Haddock, SHD; Bredeson, JV; Green, RE; Simakov, O & Rokhsar, DS
           Ancient gene linkages support ctenophores as sister to other animals. Nature (2023).
           https://doi.org/10.1038/s41586-023-05936-6

Description:
  This program is part of the Oxford Dot Plot (odp) package on github.
  This file contains a variety of functions that are useful for parsing and working with Reciprocal Best Hits (RBH) files.
  The original reason this file was written was to parse a single .rbh file into something that lists which ALGs are significantly
   located on which chromosomes.

Usage instructions:
  - See https://github.com/conchoecia/odp#getting-started
"""

import os
import pandas as pd
import sys

def hex_color_legal(hexstr) -> bool:
    """
    Checks if the hex color is legal.
    """
    if len(hexstr) != 7:
        return False
    if hexstr[0] != "#":
        return False
    for char in hexstr[1:]:
        if char not in "0123456789abcdefABCDEF":
            return False
    return True

def parse_rbh(rbhfilepath) -> pd.DataFrame:
    """
    Parses an rbh file into a pandas dataframe.
    Performs some checks to make sure that the file is legal.
      - Starts by checking that the file exists.
      - Check that there is a 'rbh' column
      - Continues to check that all the samples have a "_scaf", "_gene", and "_pos" column.
    """
    # first check that the rbhfilepath exists
    if not os.path.exists(rbhfilepath):
        raise IOError(f"The file {rbhfilepath} does not exist.")
    # now we read in the file
    df = pd.read_csv(rbhfilepath, sep = "\t")

    # Check that there is an rbh column
    if "rbh" not in df.columns:
        raise IOError(f"The rbh file, {rbhfilepath} does not have a column named 'rbh'")

    samples = [x.split("_")[0] for x in df.columns if "_scaf" in x]
    # Now check that all of the samples have a _scaf, _gene, and _pos column
    missing_fields = []
    for thissample in samples:
        for suffix in ["_scaf", "_gene", "_pos"]:
            thiscol = f"{thissample}{suffix}"
            if thiscol not in df.columns:
                missing_fields.append(thiscol)
    # print an error if there are missing fields for a sample
    if len(missing_fields) > 0:
        raise IOError(f"Error: The following columns are missing from the rbh file, {rbhfilepath} {missing_fields}. Exiting.")

    # check that all the hex colors are legal, if present in the "color" column
    if "color" in df.columns:
        for color in [x for x in df["color"] if pd.notna(x) ]:
            if not hex_color_legal(color):
                raise IOError(f"Error: The color {color} is not a legal hex color. Exiting.")
    # If we have not raised an error, then we return the dataframe
    return df

def rbh_to_scafnum(df, samplename) -> int:
    """
    Looks in an rbh file and returns the number of scaffolds in the file.
    If you know that the assembly only has chromosomes, this is a way that probably returns the chromosome number.
    """
    scafcol = f"{samplename}_scaf"
    if scafcol not in df.columns:
        raise IOError(f"The rbh file does not have a column named {scafcol}. Exiting.")
    return df[scafcol].nunique()

def rbhdf_to_alglocdf(df, minsig, ALGname) -> (pd.DataFrame, str):
    """
    This takes a .rbh filepath and returns a dataframe of the ALGs and their locations.

    Returns:
      - a tuple that contains the dataframe of the ALG coloc significance, and the sample name
      - (df, samplename)
    """
    # check that the user is inputting a dataframe. If not, they are likely inputting a filepath.
    if not isinstance(df, pd.DataFrame):
        raise IOError(f"The input to rbhdf_to_alglocdf should be a pandas dataframe. You input a {type(df)}. Did you inadvertently input a filepath?")

    # make sure that the genegroup column is present
    if "gene_group" not in df.columns:
        raise IOError(f"The rbh file, {rbhfilepath} does not have a column named 'gene_group'")

    # we need to get the sample names
    samples = [x.split("_")[0] for x in df.columns if "_scaf" in x]
    # We need to check that the ALGname is in the samples
    if not ALGname in samples:
        raise IOError(f"The ALGname, {ALGname} is not in the samples. Exiting.")
    # We need to check that the ALGname is not a part of the other sample names.
    # For example, the ALGname like BCnSSimakov2022 would be a part of a samplename BCnSSimakov2022Hydra
    othersamples = [x for x in samples if ALGname != x]
    for othersample in othersamples:
        if ALGname in othersample:
            em =  f"The ALGname, {ALGname} is a part of the samplename {othersample}."
            em += f" You shouldn't have any name overlaps. Exiting."
            raise IOError(em)
    # Right now, we are working on the assumption that we are looking specifically at the localization
    #  of a single sample with a single ALG set. There should not be more than one sample left over in
    #  othersamples.
    if len(othersamples) > 1:
        rm =  f"Error: There is more than one sample in the rbh file aside from the ALGname, {ALGname}."
        rm += f" The intended usage of this function is to use the output of the odp pipeline and"
        rm += f"  the config.yaml option: `plot_LGs: True`."
        rm += f" The rbh files you are looking for for this function will be in the directory odp/step2-figures/ALG-species_plots/ ."
        rm += f" Exiting."
        raise IOError(rm)

    # now we should get the sample name, it is the prefix to _scaf in the column that isn't the ALG's _scaf column
    samplename = othersamples[0]
    samplescafcol = f"{samplename}_scaf"
    # use chromnum as the length of the unique entries in the samplescafcol
    chromnum = len(df[samplescafcol].unique())
    # get all the rows for which whole_FET is leq than the minimum sig value
    tempdf = df[df["whole_FET"] <= minsig]
    # First we check that grouping the "gene_group" and samplescafcol should have the same whole_FET value for each row.
    # If not, this indicates that the Fisher's exact test was not being performed correctly.
    db = tempdf.groupby(["gene_group", samplescafcol])
    for name, group in db:
        if len(group["whole_FET"].unique()) != 1:
            raise IOError(f"The whole_FET column should have the same value for all rows in the same gene_group and {samplescafcol}. Exiting.")
    # Everything is fine, so we can groupby all three at the same time.
    # groupby the gene_group and the samplescafcol
    gb = tempdf.groupby(["gene_group", samplescafcol, "whole_FET"])
    # print the name of all the groups
    entries = []
    for name, group in gb:
        entries.append({"sample": samplename,
                        "gene_group": name[0],
                        "scaffold": name[1],
                        "pvalue": name[2],
                        "num_genes": len(group),
                        "frac_of_this_ALG_on_this_scaffold": len(group)/len(df[df["gene_group"] == name[0]])
                        })
    splitsdf = pd.DataFrame(entries)
    return splitsdf, samplename

def parse_ALG_rbh_to_colordf(rbh_file) -> pd.DataFrame:
    """
    Reads in an rbh file from ALGs and returns the following dataframe:

     ALGname  Color     Size
     Qb       #C72480   12
     Qc       #DCC0F3   14
     C2       #387FB2   18
     Qd       #94C47F   22
      R       #F4B93E   24
     Qa       #78A6AF   30
     A2       #8B4E67   41
     B3       #FA9A26   46
     O2       #AB5BA8   46
     Eb       #B76BED   47
    A1b       #C33D53   51
     J1       #54AB53   54
     O1       #FBD76C   55
     J2       #E64657   66
      P       #C33E51   78
     B2       #1F779A   86
      I       #3425FB   90
     B1       #2F54E3   95
      M       #A45530  102
      L       #7DC29F  104
      N       #D8BE3C  107
     Ea       #AB7E26  115
      K       #170B88  119
      H       #F04C08  135
      G       #E97B4A  138
     C1       #B07DF4  142
      F       #9B6870  145
      D       #47957F  172
    A1a       #4DB5E3  207

    The dataframe will be used later to determine plotting parameters.
    """
    df = pd.read_csv(rbh_file, sep='\t')
    # First make sure that there are columns called "gene_group", "color". We will only use these.
    if not ("gene_group" in df.columns and "color" in df.columns):
        raise IOError("The rbh file {} does not have the correct columns".format(rbh_file))
    # just subset these columns since we only need them
    df = df[["gene_group", "color"]]
    sizemap = df.groupby("gene_group").size()
    # groupby the "gene_group" column, then reconstitute with the size and the most common color for that group
    df = df.groupby("gene_group").agg({"color": lambda x: x.value_counts().index[0]})
    df = df.reset_index()
    # use the sizemap to get the size of each ALG. Use the gene_group as the key to look up the value in sizemap
    df["size"] = df["gene_group"].map(sizemap)
    # rename the columns {"gene_group": "ALGname", "color": "Color", "size": "Size"}
    df = df.rename(columns={"gene_group": "ALGname", "color": "Color", "size": "Size"})
    # sort by increasing size
    df = df.sort_values(by="Size", ascending=True)
    df = df.reset_index(drop=True)
    return df
