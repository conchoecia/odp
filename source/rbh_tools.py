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

List of functions:
    - hex_color_legal(hexstr) -> bool
    - parse_rbh(rbhfilepath) -> pd.DataFrame
    - combine_rbh(rbh_filepath1, rbh_filepath2) -> pd.DataFrame
    - rbh_to_scafnum(df, samplename) -> int
    - rbhdf_to_alglocdf(df, minsig, ALGname) -> (pd.DataFrame, str)
    - parse_ALG_rbh_to_colordf(rbh_file) -> pd.DataFrame
"""

import os
import pandas as pd
import sys
import warnings

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

    euppa-2302672-GCA964059275.1_D
    octbi-37653-GCA001194135.2_D
    Mandatory columns, dtypes, pandas column type:
      - rbh           - str    - object dtype
      - color         - str    - object dtype
      - gene_group    - str    - object dtype
      - *_gene        - str    - object dtype
      - *_scaf        - str    - object dtype
      - *_pos         - int    - int64 dtype
    Optional columns and dtypes:
      - *_breakchrom  - str    - object dtype
      - *_ix          - int    - int64 dtype
      - *_break_ix    - int    - int64 dtype
      - *_FET         - float  - float64 dtype
      - *_D           - float  - float64 dtype
    """
    # first check that the rbhfilepath exists
    if not os.path.exists(rbhfilepath):
        raise IOError(f"The file {rbhfilepath} does not exist.")

    # Now we read in the file. Anticipate that some of the columns will have NaN values.
    # If the column DOES have NaN values, then raise an error and tell the user they should just be empty.
    # Just do a try/except block to catch the error.
    # Sometimes this just breaks with large rbh (tsv) files because of how files are loaded in chunks.
    # This will need to be fixed eventually. TODO
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("error")  # Treat warnings as exceptions within this block
            df = pd.read_csv(rbhfilepath, sep="\t")
    except Warning as w:
        print(f"Warning caught during file reading. Likely because something wrong with your file {rbhfilepath}: {w}")

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
    # Make sure that the other mandatory fields are present
    other_mandatory_fields = ["gene_group", "rbh"]
    for field in other_mandatory_fields:
        if field not in df.columns:
            missing_fields.append(field)
    # print an error if there are missing fields for a sample
    if len(missing_fields) > 0:
        raise IOError(f"Error: The following mandatory columns are missing from the rbh file, {rbhfilepath} {missing_fields}. Exiting.")

    # check that all the hex colors are legal, if present in the "color" column
    if "color" in df.columns:
        for color in [x for x in df["color"] if pd.notna(x) ]:
            if not hex_color_legal(color):
                raise IOError(f"Error: The color {color} is not a legal hex color. Exiting.")

    # We need to set the dtype of different columns
    # We know a few column names that should be present in the dataframe
    object_type_columns = ["rbh", "color", "gene_group"] + \
                          [f"{x}_gene" for x in samples] + \
                          [f"{x}_scaf" for x in samples] + \
                          [f"{x}_breakchrom" for x in samples]
    for thiscol in object_type_columns:
        if thiscol in df.columns:
            df[thiscol] = df[thiscol].astype(object)
    # ints
    int_type_columns    = [f"{x}_pos" for x in samples]
    for thiscol in int_type_columns:
        if thiscol in df.columns:
            df[thiscol] = df[thiscol].astype("Int64")
    # floats
    float_type_columns  = [f"{x}_FET" for x in samples] + \
                          [f"{x}_D" for x in samples]
    for thiscol in float_type_columns:
        if thiscol in df.columns:
            df[thiscol] = df[thiscol].astype(float)

    # If we have not raised an error, then we return the dataframe
    return df

def combine_rbh_db(rbh_filepath1, rbh_filepath2) -> pd.DataFrame:
    """
    This combines two rbh database files.
    In this file type, we only have these columns:

    rbh                                             gene_group              color
    Lachesis_group13:75073921_75074550-eupbe6       Lachesis_group13        #FF7F00
    Lachesis_group22:8697662_8698153-eupbe17        Lachesis_group22        #FDBF6F
    Lachesis_group18:41128667_41129120-eupbe25      Lachesis_group18        #FFC233
    Lachesis_group21:100568739_100569082-eupbe34    Lachesis_group21        #80B1D3
    Lachesis_group1:23733546_23733930-eupbe40       Lachesis_group1         #FB9A99
    Lachesis_group15:42710712_42711059-eupbe55      Lachesis_group15        #BC80BD
    """
    # check that the two files exist
    for thisfile in [rbh_filepath1, rbh_filepath2]:
        if not os.path.exists(thisfile):
            raise IOError(f"The file {thisfile} does not exist.")

    # parse the two files
    dfs = [parse_rbh(rbh_filepath1),
           parse_rbh(rbh_filepath2)]

    # check that both dataframes have the columns "rbh", "gene_group"
    for thiscol in ["rbh", "gene_group"]:
        for thisfile in [rbh_filepath1, rbh_filepath2]:
            df = parse_rbh(thisfile)
            if thiscol not in df.columns:
                raise IOError(f"The rbh file, {thisfile} does not have a column named {thiscol}")

    for thisdf in dfs:
        # If there is no "color" column in either dataframe,
        #   then add it and make the color black
        if "color" not in thisdf.columns:
            thisdf["color"] = "#000000"
        # If there is a "color" column, then fill missing values to black
        else:
            thisdf["color"] = thisdf["color"].fillna("#000000")

    # remove all columns except "rbh", "gene_group", "color"
    for i in range(len(dfs)):
        dfs[i] = dfs[i][["rbh", "gene_group", "color"]]

    # Make sure that there are no duplicates in the two rbh columns
    for i in range(len(dfs)):
        for j in range(i+1, len(dfs)):
            if len(set(dfs[i]["rbh"]).intersection(set(dfs[j]["rbh"]))) > 0:
                raise IOError(f"The two rbh files, {rbh_filepath1} and {rbh_filepath2} have shared entries in the 'rbh' column. Exiting.")

    # merge the two dataframes
    df = pd.concat(dfs, ignore_index=True)
    return df


def combine_rbh(rbh_filepath1, rbh_filepath2) -> pd.DataFrame:
    """
    Description:
      This takes two rbh files, probably from different sources of data, and combines them into a single dataframe.

      The use case this was written for was to combine rbh files that come from different analyses,
       for example one from CNEs and one from a set of ALGs.

      Uses the columns from the first rbh file as the columns + sample for the dataframe.
       (Basically a left join)

    Column Information:
      Columns that will be removed from both dataframes:
        - *_FET         - Will be removed since these values will be outdated.
        - *_D           - Values will be outdated after merge.
        - *_ix          - No longer relevant
        - *_break_ix    - No longer relevant

      Optional columns that may be present in none, one, or both of the dataframes:
        - color         - str    - object dtype
        - *_breakchrom  - str    - object dtype

      Mandatory columns that both dataframes must have, dtypes, pandas column type:
        - rbh           - str    - object dtype
        - gene_group    - str    - object dtype
        - *_gene        - str    - object dtype
        - *_scaf        - str    - object dtype
        - *_pos         - int    - int64 dtype
    """
    # read in the two rbh files. This ensures that the files exist and are legal.
    df1 = parse_rbh(rbh_filepath1)
    df2 = parse_rbh(rbh_filepath2)

    # Make sure that the samples present in the first file are present in the second file.
    # There should only be one shared sample.
    df1_samples = sorted([x.split("_")[0] for x in df1.columns if "_scaf" in x])
    df2_samples = sorted([x.split("_")[0] for x in df2.columns if "_scaf" in x])
    # make sure there are only two samples in the two files
    if len(df1_samples) not in (1,2):
        raise IOError(f"The first rbh file, {rbh_filepath1} does not have exactly one or two samples. Exiting.")
    if len(df2_samples) not in (1,2):
        raise IOError(f"The second rbh file, {rbh_filepath2} does not have exactly one or two samples. Exiting.")
    # make sure that the two files have exactly one shared sample
    shared = set(df1_samples).intersection(set(df2_samples))
    if len(shared) != 1:
        raise IOError(f"The two rbh files, {rbh_filepath1} and {rbh_filepath2} do not have exactly one shared sample. Exiting.")
    shared_sample = shared.pop()

    # remove the columns that are not needed
    remove_cols = ["_FET", "_D", "_break_ix", "_ix",
                   "_plotindex", "_breakchrom", "_plotpos"]
    all_columns = list(set(df1.columns).union(set(df2.columns)))
    for thiscol in list(df1.columns) + list(df2.columns):
        for thissuffix in remove_cols:
            if thiscol.endswith(thissuffix):
                if thiscol in df1.columns:
                    df1 = df1.drop(thiscol, axis=1)
                if thiscol in df2.columns:
                    df2 = df2.drop(thiscol, axis=1)

    # get the unique samples in each file
    df1_unique_list = [x for x in df1_samples if x != shared_sample]
    # This list should just be one element long
    if len(df1_unique_list) != 1:
        raise IOError(f"The df1_unique list should only have one element. It has {len(df1_unique)} elements. Exiting.")
    df2_unique_list = [x for x in df2_samples if x != shared_sample]
    # This list should just be one element long
    if len(df2_unique_list) != 1:
        raise IOError(f"The df2_unique list should only have one element. It has {len(df2_unique)} elements. Exiting.")
    df1_unique = df1_unique_list[0]
    df2_unique = df2_unique_list[0]

    # Now we need to rename the columns that are unique to each file
    for thiscol in df1.columns:
        if ("_" in thiscol) and (thiscol.split("_")[-1] in ["scaf", "gene", "pos"]):
            newcolname = thiscol.replace(df1_unique, "MergedCol")
            df1 = df1.rename(columns={thiscol: newcolname})
    for thiscol in df2.columns:
        if ("_" in thiscol) and (thiscol.split("_")[-1] in ["scaf", "gene", "pos"]):
            newcolname = thiscol.replace(df2_unique, "MergedCol")
            df2 = df2.rename(columns={thiscol: newcolname})

    # Because we are merging, we need to make sure certain columns don't have shared strings.
    #  This is because something like a CNE and an ALG should not share the same label.
    #  The point of this whole function is to merge two unique forms of data into one analysis
    #  framework.
    for thiscol in ["rbh", "MergedCol_gene", f"{shared_sample}_gene"]:
        # Make sure there are no matching entries in these two columns in df1 and df2
        if len(set(df1[thiscol]).intersection(set(df2[thiscol]))) > 0:
            raise IOError(f"The columns {thiscol} in df1 and df2 have shared entries. Exiting.")

    # There are a few columns that SHOULD have matching values. Specifically, the columns
    #  with the scaffolds of the shared sample should have matching values, meaning that these
    #  markers fall on the same scaffolds. This is because the point of doing this is that
    #  merging the rbh files is to get better resolution of markers.
    for thiscol in [f"{shared_sample}_scaf"]:
        # Make sure there is at least one matching entry. That is exceptionally forgiving
        #  but at least it is a mark that these are derived from the same genome assembly.
        if len(set(df1[thiscol]).intersection(set(df2[thiscol]))) == 0:
            raise IOError(f"The columns {thiscol} in df1 and df2 have no shared entries. Exiting.")

    # Now we can merge the two dataframes. We just stack them on top of each other, then sort by
    #  the shared sample scaffold then position
    #print("df1 before merge: \n", df1)
    #print("df2 before merge: \n", df2)
    df = pd.concat([df1, df2], ignore_index=True)
    df = df.sort_values(by=[f"{shared_sample}_scaf", f"{shared_sample}_pos"]).reset_index(drop=True)
    #print()
    #print("df after merge: \n", df)
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
