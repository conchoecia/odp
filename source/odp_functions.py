"""
These are functions that are shared by odp, odp_trio, and other scripts
"""

# this is all needed to load our custom fasta parser
import gzip
import os
import sys
snakefile_path = os.path.dirname(os.path.abspath(__file__))
dependencies_path = os.path.join(snakefile_path, "../dependencies")
sys.path.insert(1, dependencies_path)
import afp as fasta

# ODP-specific imports
import odp_plotting_functions as odp_plot

# other standard python libraries
from itertools import combinations
from itertools import groupby
from itertools import product
from operator import itemgetter
import hashlib

# non-standard dependencies
import pandas as pd
import numpy as np

# TODO implement this function soon after release
#def gff_or_chrom(filepath):
#    """
#    This function takes a filepath and determines if it is a gff or a chrom file.
#    Returns "gff" or "chrom", or throws an error if the software can't decide.
#    """
#    # First check what the file extension is

def hmmsearch_get_mem_mb(wildcards, attempt):
    """
    The amount of RAM needed for miniprot is highly variable.
    """
    attemptdict = {1: 2000,
                   2: 4000,
                   3: 8000,
                   4: 16000,
                   5: 32000,
                   6: 64000,
                   7: 128000,
                   8: 256000}
    return attemptdict[attempt]


def filthmm_get_mem_mb(wildcards, attempt):
    """
    The amount of RAM needed for miniprot is highly variable.
    """
    attemptdict = {1: 1000,
                   2: 2000,
                   3: 4000,
                   4: 8000,
                   5: 16000,
                   6: 320000,
                   7: 640000}
    return attemptdict[attempt]


def gzhmm_get_mem_mb(wildcards, attempt):
    """
    The amount of RAM needed for miniprot is highly variable.
    """
    attemptdict = {1: 1000,
                   2: 2000,
                   3: 4000,
                   4: 8000,
                   5: 16000,
                   6: 320000,
                   7: 640000}
    return attemptdict.get(attempt, attemptdict[max(attemptdict)])



def hmm_against_prots_get_mem_mb(wildcards, attempt):
    """
    The amount of RAM needed for filtering prots could change.
    """
    attemptdict = {1: 5000,
                   2: 16000,
                   3: 64000,
                   4: 256000}
    return attemptdict[attempt]

def tmp_unzip_get_mem_mb(wildcards, attempt):
    """
    The amount of RAM needed for miniprot is highly variable.
    """
    attemptdict = {1: 2000,
                   2: 4000,
                   3: 8000,
                   4: 16000}
    return attemptdict[attempt]


def chromsize_to_s2c2s(chromsize_path):
    """
    this reads in a chromsize file and returns a dictionary of species to scaffold to scaflen
    """
    s2c2s = {}
    with open(chromsize_path, "r") as f:
        for line in f:
            line = line.strip()
            if line:
                thissp, thisscaf, thislen = line.split("\t")
                if thissp not in s2c2s:
                    s2c2s[thissp] = {}
                s2c2s[thissp][thisscaf] = int(thislen)
    return s2c2s

def calc_D_for_y_and_x(
    df,
    xsample=None,
    ysample=None,
    x_offset=None,
    y_offset=None,
    x_scaf_to_len=None,
    y_scaf_to_len=None,
    **kwargs,
):
    """Calculate D for both the x and y axes.

    This function supports two different dataframe structures:

    1. Break-based mode where ``xsample`` and ``ysample`` are provided and the
       dataframe contains ``*_breakchrom`` and ``*_break_ix`` columns.
    2. Scaffold-based mode where offset dictionaries and scaffold length
       mappings are supplied for both axes.
    """

    # Break-based mode
    if (
        xsample is not None
        and ysample is not None
        and f"{xsample}_breakchrom" in df.columns
        and f"{ysample}_breakchrom" in df.columns
    ):
        for thisdir, oppositexy in [(xsample, ysample), (ysample, xsample)]:
            df = df.sort_values(
                by=[f"{thisdir}_breakchrom", f"{thisdir}_break_ix"], ascending=True
            ).reset_index(drop=True)
            breaks = df[f"{thisdir}_breakchrom"].unique()
            thisdir_dfs = []
            for thisx in breaks:
                xdf = (
                    df.loc[df[f"{thisdir}_breakchrom"] == thisx,].copy().reset_index(drop=True)
                )
                df2 = pd.get_dummies(xdf[f"{oppositexy}_breakchrom"])
                df2_xiL = df2.apply(lambda x: x.rolling(20).mean(), axis=0)
                df2_xiR = (
                    df2.apply(lambda x: x.iloc[::-1].rolling(20).mean(), axis=0)
                    .iloc[::-1]
                    .set_index(df2.index - 1)
                    .iloc[1:]
                )
                subtractdf = df2_xiR.fillna(0) - df2_xiL.fillna(0)
                xdf[f"{thisdir}_D"] = subtractdf.apply(
                    lambda x: np.sqrt(np.square(x).sum()), axis=1
                )
                thisdir_dfs.append(xdf)
            df = pd.concat(thisdir_dfs)
        df.reset_index(drop=True, inplace=True)
        return df

    # Scaffold-based mode
    elif (
        x_offset is not None
        and y_offset is not None
        and x_scaf_to_len is not None
        and y_scaf_to_len is not None
    ):
        df = df.dropna()
        for thisdir in ["x", "y"]:
            df = df.sort_values(by=[f"{thisdir}middle"])
            df.reset_index(drop=True, inplace=True)

            unique_x = df[f"{thisdir}scaf"].unique()
            thisdir_dfs = []
            for thisx in unique_x:
                xdf = df.loc[df[f"{thisdir}scaf"] == thisx,].copy().reset_index(drop=True)
                oppositexy = "y" if thisdir == "x" else "x"
                this_offset = x_offset if thisdir == "x" else y_offset
                this_scaf_to_len = (
                    x_scaf_to_len if thisdir == "x" else y_scaf_to_len
                )
                df2 = pd.get_dummies(xdf[f"{oppositexy}scaf"])
                df2_xiL = df2.apply(lambda x: x.rolling(20).mean(), axis=0)
                df2_xiR = (
                    df2.apply(lambda x: x.iloc[::-1].rolling(20).mean(), axis=0)
                    .iloc[::-1]
                    .set_index(df2.index - 1)
                    .iloc[1:]
                )
                subtractdf = df2_xiR.fillna(0) - df2_xiL.fillna(0)
                D = subtractdf.apply(lambda x: np.sqrt(np.square(x).sum()), axis=1)
                xdf[f"D{thisdir}"] = D
                # Float dtype: barmiddle is start + width/2 which can be
                # non-integer; pandas 3.x makes the implicit cast an error.
                xdf[f"D{thisdir}_barleft"] = 0.0
                xdf[f"D{thisdir}_barmiddle"] = 0.0
                xdf[f"D{thisdir}_barright"] = 0.0
                xdf[f"D{thisdir}_barwidth"] = 0.0
                for i, row in xdf.iterrows():
                    barleft = -1
                    barright = -1
                    if len(xdf) > 1:
                        if i == 0:
                            thisend = row[f"{thisdir}stop"]
                            nextstart = xdf.loc[i + 1, f"{thisdir}start"]
                            barleft = this_offset[thisx]
                            barright = thisend + ((nextstart - thisend) / 2)
                        elif i == (len(xdf) - 1):
                            prevend = xdf.loc[i - 1, f"{thisdir}stop"]
                            thisstart = row[f"{thisdir}start"]
                            barleft = prevend + ((thisstart - prevend) / 2)
                            barright = this_scaf_to_len[thisx]
                        else:
                            prevend = xdf.loc[i - 1, f"{thisdir}stop"]
                            thisstart = row[f"{thisdir}start"]
                            thisend = row[f"{thisdir}stop"]
                            nextstart = xdf.loc[i + 1, f"{thisdir}start"]
                            barleft = prevend + ((thisstart - prevend) / 2)
                            barright = thisend + ((nextstart - thisend) / 2)
                    xdf.loc[i, f"D{thisdir}_barleft"] = barleft
                    xdf.loc[i, f"D{thisdir}_barright"] = barright
                    xdf.loc[i, f"D{thisdir}_barmiddle"] = barleft + ((barright - barleft) / 2)
                    xdf.loc[i, f"D{thisdir}_barwidth"] = barright - barleft + 1
                thisdir_dfs.append(xdf)
            df = pd.concat(thisdir_dfs)
        df = df.sort_values(by=["xmiddle"])
        df.reset_index(drop=True, inplace=True)
        return df

    else:
        raise ValueError("Insufficient arguments to calculate D for x and y")


def filtprot_get_mem_mb(wildcards, attempt):
    """
    The amount of RAM needed for filtering prots could change.
    """
    attemptdict = {1: 2000,
                   2: 4000,
                   3: 16000,
                   4: 64000,
                   5: 256000,
                   6: 512000}
    return attemptdict.get(attempt, max(attemptdict.values()))

def general_legal_run():
    """
    imports required:
      - os
      - sys

    Checks if the run itself is legal. We need to check for:

    1. This program is not being run in a subdirectory of the odp install.
       We do not allow this, as some of the outfiles may overwrite program files.
    """
    snakefile_path = os.path.dirname(os.path.abspath(__file__))
    odp_path       = os.path.abspath(os.path.join(snakefile_path, ".."))
    safe_dirs      = [os.path.join(odp_path, "tests/test_odp_basic")]
    cwd            = os.path.abspath(os.getcwd())

    # if we are in the odp directory, but in the test directory, that's fine
    # test if we are in the odp directory
    crash = False
    if odp_path in cwd:
        crash = True
    if cwd in safe_dirs:
        # we don't care, actually. Just put the flag back to False
        print(f"We are in the test directory, {cwd}", file = sys.stderr)
        crash = False

    # check to see if we crash
    if crash:
        # raise an error telling the user not to run the analysis in the odp directory
        outmessage =  "*********************************************************************\n"
        outmessage += "* ERROR:\n"
        outmessage += "*  You are running this program in the odp install directory.\n"
        outmessage += "*  The directory where odp is installed is: " + odp_path + "\n"
        outmessage += "*  The directory where this analysis is being run is: " + cwd + "\n"
        outmessage += "*\n"
        outmessage += "*  The reason this is problematic is that some of the output files\n"
        outmessage += "*   may overwrite program files.\n"
        outmessage += "*\n"
        outmessage += "*  Please run this analysis in a different directory.\n"
        outmessage += "*********************************************************************\n"
        # now use this message for the error and exit the program
        raise ValueError(outmessage)

def reciprocal_best_permissive_blastp_or_diamond_blastp(
        x_to_y_blastp_results, y_to_x_blastp_results, outfile):
    """
    This function finds reciprocal best blastp hits between two samples.
    The input is a blastp results file where x was blasted against y,
      and a blastp results file where y was blasted against x.

    The output format is just the rows of the blastp results from the x_to_y file.
    Saves it as a df to outfile.

    This algorithm does not have an absolute best, but leaves all possible
      best hits based on e-value to be filtered out later by
      analyzing a graph of the blast results
    """
    f_raw = pd.read_csv(x_to_y_blastp_results, sep="\t")
    f_raw.columns = ["qseqid", "sseqid", "pident", "length",
                   "mismatch", "gapopen", "qstart", "qend",
                   "sstart", "send", "evalue", "bitscore"]
    # Keep every row whose evalue equals the per-qseqid minimum.
    # Vectorised transform avoids the groupby-apply pattern that pandas
    # 3.x drops grouping columns from.
    fdf = f_raw.loc[
        f_raw.groupby("qseqid")["evalue"].transform("min") == f_raw["evalue"]
    ].reset_index(drop=True)

    r_raw = pd.read_csv(y_to_x_blastp_results, sep="\t")
    r_raw.columns = ["qseqid", "sseqid", "pident", "length",
                   "mismatch", "gapopen", "qstart", "qend",
                   "sstart", "send", "evalue", "bitscore"]
    rdf = r_raw.loc[
        r_raw.groupby("qseqid")["evalue"].transform("min") == r_raw["evalue"]
    ].reset_index(drop=True)
    rdf.columns = ["sseqid", "qseqid", "pident", "length",
                   "mismatch", "gapopen", "qstart", "qend",
                   "sstart", "send", "evalue", "bitscore"]
    rdf = rdf[["qseqid","sseqid"]]

    #These are the singleton RBH
    new_df = pd.merge(fdf, rdf,  how='inner', left_on=['qseqid','sseqid'], right_on = ['qseqid','sseqid'])
    new_df.to_csv(outfile, sep="\t", index = False, header = False)

def reciprocal_best_hits_blastp_or_diamond_blastp(
        x_to_y_blastp_results, y_to_x_blastp_results, outfile):
    """
    This function finds reciprocal best blastp hits between two samples.
    The input is a blastp results file where x was blasted against y,
      and a blastp results file where y was blasted against x.

    The output format is just the rows of the blastp results from the x_to_y file.
    Saves it as a df to outfile.

    This algorithm is permissive in that it finds the best hits between the two
      species even if the e-values for the "best hit" are equivalent. This fixes
      one of the problems with blastp results. The results are still reciprocal
      best, though.
    """
    f_raw = pd.read_csv(x_to_y_blastp_results, sep="\t")
    f_raw.columns = ["qseqid", "sseqid", "pident", "length",
                   "mismatch", "gapopen", "qstart", "qend",
                   "sstart", "send", "evalue", "bitscore"]
    fdf = f_raw.sort_values(["qseqid", "bitscore", "evalue", "pident"], ascending=[True, False, True, False]).drop_duplicates(subset="qseqid")

    r_raw = pd.read_csv(y_to_x_blastp_results, sep="\t")
    r_raw.columns = ["qseqid", "sseqid", "pident", "length",
                   "mismatch", "gapopen", "qstart", "qend",
                   "sstart", "send", "evalue", "bitscore"]
    rdf = r_raw.sort_values(["qseqid", "bitscore", "evalue", "pident"], ascending=[True, False, True, False]).drop_duplicates(subset="qseqid")
    rdf.columns = ["sseqid", "qseqid", "pident", "length",
                   "mismatch", "gapopen", "qstart", "qend",
                   "sstart", "send", "evalue", "bitscore"]
    rdf = rdf[["sseqid","qseqid"]]

    #These are the singleton RBH
    new_df = pd.merge(fdf, rdf,  how='inner', left_on=['qseqid','sseqid'], right_on = ['qseqid','sseqid'])
    #these rows are a little pedantic and we don't really need to do them
    new_df = new_df.sort_values(["qseqid","bitscore"],
                                ascending=[True, False]).drop_duplicates(
                                    subset="qseqid")
    new_df = new_df.sort_values(["sseqid","bitscore"],
                                ascending=[True, False]).drop_duplicates(
                                    subset="sseqid")

    # now filter
    f_seqs = new_df["qseqid"]
    r_seqs = new_df["sseqid"]
    # Keep every row whose evalue equals the per-qseqid minimum. The
    # vectorised `transform("min")` form is pandas-3 compatible (no
    # groupby-apply deprecation) and avoids the column-drop ambiguity.
    fdf = f_raw.copy()
    fdf = fdf.loc[
        fdf.groupby("qseqid")["evalue"].transform("min") == fdf["evalue"]
    ].reset_index(drop=True)
    rdf = r_raw.copy()
    rdf = rdf.loc[
        rdf.groupby("qseqid")["evalue"].transform("min") == rdf["evalue"]
    ].reset_index(drop=True)

    # only get the things that we haven't seen yet
    fdf2 = fdf.loc[~fdf["qseqid"].isin(f_seqs)]
    fdf2 = fdf2.loc[~fdf2["sseqid"].isin(r_seqs)]
    rdf2 = rdf.loc[ ~rdf["sseqid"].isin(f_seqs)]
    rdf2 = rdf2.loc[~rdf2["qseqid"].isin(r_seqs)]

    #swap columns for merge
    rdf2.columns = ["sseqid", "qseqid", "pident", "length",
                   "mismatch", "gapopen", "qstart", "qend",
                   "sstart", "send", "evalue", "bitscore"]
    rdf2 = rdf2[["sseqid","qseqid"]]

    new_df2 = pd.merge(fdf2, rdf2,  how='inner', left_on=['qseqid','sseqid'], right_on = ['qseqid','sseqid'])
    # get rid of duplicates
    new_df2 = new_df2.sort_values(["qseqid","bitscore"],
                                ascending=[True, False]).drop_duplicates(
                                    subset="qseqid")
    new_df2 = new_df2.sort_values(["sseqid","bitscore"],
                                ascending=[True, False]).drop_duplicates(
                                    subset="sseqid")

    # this is also pedantic and shouldn't do anything
    finaldf = pd.concat([new_df, new_df2])
    prelen = len(finaldf)
    finaldf = finaldf.sort_values(["qseqid","bitscore"],
                                ascending=[True, False]).drop_duplicates(
                                    subset="qseqid")
    finaldf = finaldf.sort_values(["sseqid","bitscore"],
                                ascending=[True, False]).drop_duplicates(
                                    subset="sseqid")
    if prelen != len(finaldf):
        raise IOError("something happened in parsing that shouldn't have. These filtering steps should not have done anything")
    finaldf.to_csv(outfile, sep="\t", index = False, header = False)

def check_file_exists(filepath) -> bool:
    """
    checks if a file exists.
    If not, raises an error
    """
    if not os.path.isfile(filepath):
        outmessage =  "*********************************************************************\n"
        outmessage += "* ERROR:\n"
        outmessage += "*  This file does not exist: " + str(filepath) + "\n"
        outmessage += "*********************************************************************\n"
        raise IOError(outmessage)
    else:
        return True

def chrom_file_is_legal(chrompath):
    """
    Checks if a chrom file is legal. Can be gzipped or not.
    Columns (and types) are:
     - protein_id (string)
     - scaffold (string)
     - strand (string)
     - start (int)
     - stop (int)

    This doesn't check if the proteins or scaffolds are legal. It simply checks if the file format is legal.

    The function that checks if this matches the protein and genome file is check_species_input_legality().

    BREAK CONDITIONS
      - Any of the fields have leading or trailing whitespace
      - Any of the following strings appear in theses respective columns: pid     scaf    strand  start   stop
      - Field 2 isn't a ['+', '-', '.']
      - Field 3 can't be converted to an int
      - Field 4 can't be converted to an int

    If any of the columns don't match this, or if there is a header string, returns False.
    If everything is good, returns True.
    """
    # 1. check that the file exists
    check_file_exists(chrompath)
    # 2. Open the file for whatever type it is. afp.get_open_func returns a
    # callable (e.g. gzip.open / open); call it to get the actual handle.
    opener = fasta.get_open_func(chrompath)
    with opener(chrompath, "rt") as chromhandle:
        # go through the file line by line and inspect each element
        for line in chromhandle:
            line = line.strip()
            if line:
                fields = line.strip().split("\t")
                # check if any of the fields have leading or trailing whitespace
                for field in fields:
                    if field != field.strip():
                        print("There is leading or trailing whitespace in this field: " + field)
                        return False
                # check if any of the fields are the header strings
                if fields[0] == "pid" or fields[1] == "scaf" or fields[2] == "strand" or fields[3] == "start" or fields[4] == "stop":
                    print("One of the fields is a header string: " + str(fields))
                    return False
                # Check if field 2 is a ['+', '-', '.']
                if fields[2] not in ['+', '-', '.']:
                    print("Field 2 is not a ['+', '-', '.']: " + str(fields))
                    return False
                # check that field 3 is able to be converted to an int
                if not fields[3].isdigit():
                    print("Field 3 is not an int: " + str(fields))
                    return False
                # check that field 4 is able to be converted to an int
                if not fields[4].isdigit():
                    print("Field 4 is not an int: " + str(fields))
                    return False
    # if we get here, everything is good
    return True

def check_species_input_legality(fastapath, peppath, chrompath, dup_proteins_allowed = False) -> bool:
    """
    This function checks that the input files are legal.
    There are certain fields that are required,
      and they must be in a specific format.

    First read in the genome assembly fasta file:
      1. Check that the file exists
      2. Check that each sequence ID exists only once

    Then read in the protein file:
      1. Check that the file exists
      2. Check that each sequence ID exists only once
      3. Check that there are no duplicate protein sequences

    Lastly, read in the .chrom file:
      1. Check that the file exists
      2. Check that the proteins in column 1 were seen in the protein fasta file
      3. Check that the scaffolds were seen in the genome assembly fasta file

    """
    # we need to check that dup_proteins_allowed is type bool
    if not isinstance(dup_proteins_allowed, bool):
        raise IOError("dup_proteins_allowed must be a boolean, True or False. However, it is " + str(dup_proteins_allowed) + " of type " + str(type(dup_proteins_allowed)) + ".")

    # PARSE AND CHECK THE GENOME ASSEMBLY
    # 1. check that the file exists
    check_file_exists(fastapath)
    # 2. check that each sequence ID exists only once
    genome_headers = set()
    duplicates     = set()
    for record in fasta.parse(fastapath):
        if record.id not in genome_headers:
            genome_headers.add(record.id)
        else:
            duplicates.add(record.id)
    if len(duplicates) > 0:
        dupstring = "".join(["*    - " + str(x) + "\n" for x in sorted(duplicates)[:3]])
        # raise an error because each ID should only occur once
        outmessage =  "*********************************************************************\n"
        outmessage += "* ERROR:\n"
        outmessage += "*  There is a genome assembly with duplicate sequence headers.\n"
        outmessage += "*  Each sequence in the genome assembly must have a unique ID.\n"
        outmessage += "*\n"
        outmessage += "*  The assembly with the problem is: " + fastapath + "\n"
        outmessage += "*  There are " + str(len(duplicates)) + " duplicate sequence headers.\n"
        outmessage += "*  Here are the first 1 to 3:\n"
        outmessage += dupstring
        outmessage += "*\n"
        outmessage += "*  The reason this is problematic is that we cannot distinguish\n"
        outmessage += "*   between two separate sequences with the same header.\n"
        outmessage += "*\n"
        outmessage += "*  Please remove the duplicate sequence headers from the fasta file,\n"
        outmessage += "*   regenerate the protein fasta and chrom files, and try again.\n"
        outmessage += "*********************************************************************\n"
        raise IOError(outmessage)

    # PARSE AND CHECK THE PROTEIN FILE
    # 1. check that the file exists
    check_file_exists(peppath)
    # 2. check that each sequence ID exists only once
    protein_headers      = set()
    duplicate_headers    = set()
    sequence_hashes      = set()
    duplicate_sequences  = set()
    for record in fasta.parse(peppath):
        if record.id not in protein_headers:
            protein_headers.add(record.id)
        else:
            duplicate_headers.add(record.id)
        seq_hash = hashlib.sha256(str(record.seq).encode()).digest()
        if seq_hash not in sequence_hashes:
            sequence_hashes.add(seq_hash)
        else:
            duplicate_sequences.add(record.id)

    if len(duplicate_headers) > 0:
        dupstring = "".join(["*    - " + str(x) + "\n" for x in sorted(duplicate_headers)[:3]])
        # raise an error because each ID should only occur once
        outmessage =  "*********************************************************************\n"
        outmessage += "* ERROR:\n"
        outmessage += "*  There is a protein fasta with duplicate sequence headers.\n"
        outmessage += "*  Each sequence in the protein fasta must have a unique ID.\n"
        outmessage += "*  Otherwise, how are we to know which sequence is which?\n"
        outmessage += "*\n"
        outmessage += "*  The protein pep with the problem is: " + peppath + "\n"
        outmessage += "*  There are " + str(len(duplicate_headers)) + " duplicate sequence headers.\n"
        outmessage += "*  Here are the first 1 to 3:\n"
        outmessage += dupstring
        outmessage += "*\n"
        outmessage += "*  The reason this is problematic is that we cannot distinguish\n"
        outmessage += "*   between two separate sequences with the same header.\n"
        outmessage += "*\n"
        outmessage += "*  Please remove the duplicate sequence headers from the protein fasta\n"
        outmessage += "*   file, regenerate the chrom files, and try again.\n"
        outmessage += "*********************************************************************\n"
        raise IOError(outmessage)

    # 3. Check that each sequence ID exists only once.
    #    This is the only case where we can choose to allow duplicates.
    if not dup_proteins_allowed:
        if len(duplicate_sequences) > 0:
            dupstring = "".join(["*    - " + str(x) + "\n" for x in sorted(duplicate_sequences)[:3]])
            # raise an error because each ID should only occur once
            outmessage =  "*********************************************************************\n"
            outmessage += "* ERROR:\n"
            outmessage += "*  Some protein sequences in your file are identical.\n"
            outmessage += "*  Each protein sequence must be unique.\n"
            outmessage += "*  If you aren't sure why, PLEASE READ THESE TWO LINKS below:\n"
            outmessage += "*    - https://github.com/conchoecia/odp/issues/49\n"
            outmessage += "*    - https://github.com/conchoecia/odp/issues/62\n"
            outmessage += "*\n"
            outmessage += "*  The protein fasta with the problem is: " + peppath + "\n"
            outmessage += "*  There are " + str(len(duplicate_sequences)) + " duplicate sequences.\n"
            outmessage += "*  Here are the first 1 to 3:\n"
            outmessage += dupstring
            outmessage += "*\n"
            outmessage += "*  The reason this is problematic is that duplicate protein seqs\n"
            outmessage += "*   may interfere with proper reciprocal blastp match detection.\n"
            outmessage += "*\n"
            outmessage += "*  Please remove the identical sequences from the protein fasta\n"
            outmessage += "*   file, regenerate the chrom files, and try again.\n"
            outmessage += "*                          -- OR --\n"
            outmessage += "*  !!! IF YOU WANT THIS ERROR TO GO AWAY without modifying your data,\n"
            outmessage += "*   set the 'duplicate_proteins' line in your 'config.yaml'\n"
            outmessage += "*   to \"pass\". The default is \"fail\"\n"
            outmessage += "*  !!! In other words, add this line to your 'config.yaml' file:\n"
            outmessage += "*   ```\n"
            outmessage += "*   duplicate_proteins: \"pass\"\n"
            outmessage += "*   ```\n"
            outmessage += "*\n"
            outmessage += "*  Your final file would look something like this:\n"
            outmessage += "*   ```\n"
            outmessage += "*   ignore_autobreaks: True\n"
            outmessage += "*   diamond_or_blastp: \"diamond\"\n"
            outmessage += "*   duplicate_proteins: \"pass\"\n"
            outmessage += "*   plot_LGs: True\n"
            outmessage += "*   plot_sp_sp: True\n"
            outmessage += "*   species:\n"
            outmessage += "*     Celegans:\n"
            outmessage += "*       proteins: /path/to/proteins_in_Cel_genome.fasta\n"
            outmessage += "*       chrom: /path/to/Cel_genome_annotation.chrom\n"
            outmessage += "*       genome: /path/to/Cel_genome_assembly.fasta\n"
            outmessage += "*       minscafsize: 1000000  # Only plots scaffolds that are 1 Mbp or longer\n"
            outmessage += "*     Homosapiens:\n"
            outmessage += "*       proteins: /path/to/Human_prots.fasta\n"
            outmessage += "*       chrom: /path/to/Human_annotation.chrom\n"
            outmessage += "*       genome: /path/to/Human_genome_assembly.fasta\n"
            outmessage += "*       minscafsize: 8000000  # Only plots scaffolds that are 8 Mbp or larger\n"
            outmessage += "*   ```\n"
            outmessage += "*********************************************************************\n"
            raise IOError(outmessage)

    # PARSE AND CHECK THE CHROM FILE
    # 1. check that the file exists
    check_file_exists(chrompath)
    proteins_not_in_pep    = set()
    scaffolds_not_in_fasta = set()
    opener = fasta.get_open_func(chrompath)
    with opener(chrompath, "rt") as chromhandle:
        for line in chromhandle:
            line = line.strip()
            if line:
                fields = line.split("\t")
                # check that the protein was seen in the protein fasta file
                protid = fields[0]
                scaffold = fields[1]
                if protid not in protein_headers:
                    proteins_not_in_pep.add(protid)
                if scaffold not in genome_headers:
                    scaffolds_not_in_fasta.add(scaffold)

    # 2. Check that the proteins in column 1 were seen in the protein fasta file
    if len(proteins_not_in_pep) > 0:
        #   if there are any duplicates, print out a string of the first 3
        proteins_not_in_pep   = sorted(list(proteins_not_in_pep))
        # raise an error because the proteins should have been seen already
        dupstring = "".join(["*    - " + str(x) + "\n" for x in proteins_not_in_pep[:3]])
        outmessage =  "*********************************************************************\n"
        outmessage += "* ERROR:\n"
        outmessage += "*  Some proteins in the .chrom file were not seen in the protein\n"
        outmessage += "*   .fasta file. This is problematic, because it could indicate\n"
        outmessage += "*   missing data, a problem generating the chrom file, or a problem\n"
        outmessage += "*   generating the protein fasta file.\n"
        outmessage += "*\n"
        outmessage += "*  The chrom file with the problem is: " + chrompath + "\n"
        outmessage += "*  There are " + str(len(proteins_not_in_pep)) + " proteins in the\n"
        outmessage += "*    .chrom not seen in the protein .fasta\n"
        outmessage += "*  Here are the first 1 to 3:\n"
        outmessage += dupstring
        outmessage += "*\n"
        outmessage += "*  The reason this is problematic is that we need to access every\n"
        outmessage += "*   protein specified in the .chrom file, but it is unavailable.\n"
        outmessage += "*\n"
        outmessage += "*  For example, in you sample one of the proteins we found in the\n"
        outmessage += "*    above chrom file was: " + proteins_not_in_pep[0] + "\n"
        outmessage += "*  This means that in the protein fasta file, there needs to be an\n"
        outmessage += "*    entry with the same name.\n"
        outmessage += "*  In the protein fasta file, there should be one protein that looks\n"
        outmessage += "*    like this, with a > character, the protein ID from the chrom file,\n"
        outmessage += "*    and then a newline or a space character ' '\n"
        outmessage += "*    ```\n"
        outmessage += "*    >" + proteins_not_in_pep[0] + " Optional sequence description here\n"
        outmessage += "*    MSNKKRN... (the protein's sequences on this and subsequent lines.)\n"
        outmessage += "*    >Next_protein_sequence\n"
        outmessage += "*    MNELSKENNIE....\n"
        outmessage += "*    ```\n"
        outmessage += "*  Please investigate whether there are too many entries in the .chrom\n"
        outmessage += "*   file, or if something is missing from the protein .fasta file.\n"
        outmessage += "*   Then, fix your files and re-run this pipeline.\n"
        outmessage += "*\n"
        outmessage += "* If you need more help, please visit:\n"
        outmessage += "*  https://github.com/conchoecia/odp?tab=readme-ov-file#chrom-file-specifications\n"
        outmessage += "*\n"
        outmessage += "*********************************************************************\n"
        raise IOError(outmessage)

    # 3. Check that the scaffolds were seen in the genome assembly fasta file.
    if len(scaffolds_not_in_fasta):
        # Error. Scaffolds specified in .chrom file but missing in the genome assembly fasta.
        dupstring = "".join(["*    - " + str(x) + "\n" for x in sorted(scaffolds_not_in_fasta)[:3]])
        outmessage =  "*********************************************************************\n"
        outmessage += "* ERROR:\n"
        outmessage += "*  Some scaffolds in the .chrom file were not seen in the genome\n"
        outmessage += "*   assembly .fasta file.\n"
        outmessage += "*\n"
        outmessage += "*  The chrom file with the problem is: " + chrompath + "\n"
        outmessage += "*  There are " + str(len(scaffolds_not_in_fasta)) + " scaffolds in the .chrom not seen in the genome .fasta\n"
        outmessage += "*  Here are the first 1 to 3:\n"
        outmessage += dupstring
        outmessage += "*\n"
        outmessage += "*  The reason this is problematic is that we need to access every\n"
        outmessage += "*   scaffold specified in the .chrom file, but it is unavailable.\n"
        outmessage += "*\n"
        outmessage += "*  Please investigate whether there are too many entries in the .chrom\n"
        outmessage += "*   file, or if something is missing from the genome .fasta file.\n"
        outmessage += "*   Then, fix your files and re-run this pipeline.\n"
        outmessage += "*********************************************************************\n"
        # The error message was previously built but never raised, so
        # invalid configs slipped past undetected. Raise it now.
        raise IOError(outmessage)

    # everything passed
    return True

def check_legality(config):
    """
    This function checks for legal config entries.
    This is useful for finding misspelled entries.
    Just checks to see if the arguments in this config file are legal.
    """
    # The following strings are illegal and may have been used in previous versions of the program
    #  - "prot_to_loc"
    #  - "prot_to_loc"
    legal = ["assembly_accession", "proteins", "chrom", "genome", "genus",
             "minscafsize", "manual_breaks", "chrom_to_color",
             "plotorder", "prot_to_group", "species", "taxid" ]
    illegal = set()

    for this_axis in ["species"]:
        if this_axis in config:
            for this_sample in config[this_axis]:
                for key in config[this_axis][this_sample]:
                    if key not in legal:
                        illegal.add(key)
    if len(illegal) > 0:
        print("We found some fields in your config file that are not used by this program.")
        print("The only fields allowed for individual samples are:")
        for key in legal:
            print("  - {}".format(key))
        print("The keys that we found that are not allowed/in the list above are:")
        for key in illegal:
            print("  - {}".format(key))
        sys.exit()
    if "species" in config:
        for thissample in config["species"]:
            if "_" in thissample:
                raise IOError("Sample names can't have '_' char: {}".format(thissample))

def flatten(list_of_lists):
    """flatten a list of lists, unique only"""
    return list(set([item for sublist in list_of_lists for item in sublist]))

def expand_avoid_matching_x_and_y(filestring, xsamples, ysamples):
    """
    this works like snakemake's expand function but does not generate
     files where xsample equals ysample

    outputs a list of files
    """
    outlist = []
    for xsamp in xsamples:
        for ysamp in ysamples:
            if xsamp != ysamp:
                outlist.append(filestring.format(xsamp, ysamp))
    return outlist

def expand_avoid_matching_x_and_y_third(filestring, xsamples, ysamples, third):
    """
    this works like snakemake's expand function but does not generate
     files where xsample equals ysample

    outputs a list of files.

    Note March 20th, 2021 - the "third" aspect is related to coloring by another species
    """
    outlist = []
    for xsamp in xsamples:
        for ysamp in ysamples:
            if xsamp != ysamp:
                for t in third:
                    print("hey", xsamp, ysamp, third)
                    outlist.append(filestring.format(xsamp, ysamp, t))
    return outlist

def expand_avoid_matching(filestring, **kwargs):
    """
    this works like snakemake's expand function but does not generate
     files where xsample equals ysample

    outputs a list of files
    """
    outlist = []
    keys = [x for x in kwargs]
    values = [[y for y in kwargs[x]] for x in keys]
    nonmatching_products = [x for x in list(product(*values)) if len(set(x)) == len(x)]
    for entry in nonmatching_products:
        these_kwargs = {}
        for i in range(len(entry)):
            these_kwargs[keys[i]] = entry[i]
        outlist.append(filestring.format(**these_kwargs))
    return [x for x in outlist]

def expand_avoid_matching_all_third(filestring, **kwargs):
    """
    this works like snakemake's expand function but does not generate
     files where xsample equals ysample

    outputs a list of files
    """
    if "third" not in kwargs:
        raise IOError("third not in kwargs")
    outlist = []
    keys = [x for x in kwargs if x != "third"]
    values = [[y for y in kwargs[x]] for x in keys]
    nonmatching_products = [x for x in list(product(*values)) if len(set(x)) == len(x)]
    for entry in nonmatching_products:
        these_kwargs = {}
        for i in range(len(entry)):
            these_kwargs[keys[i]] = entry[i]
        for t in kwargs["third"]:
            temp_kwargs = these_kwargs
            temp_kwargs["third"] = t
            outlist.append(filestring.format(**temp_kwargs))
    return [x for x in outlist]


def open_text_maybe_gzip(path, encoding="utf-8"):
    with open(path, "rb") as fh:
        head = fh.read(2)
    if head == b"\x1f\x8b":
        return gzip.open(path, "rt", encoding=encoding)
    return open(path, "rt", encoding=encoding)

def filter_fasta_chrom(chrom_file, input_fasta, output_fasta):
    """
    Keep only proteins whose IDs appear in chrom_file; write to output_fasta (gz or plain).
    RAM: O(#unique IDs in chrom) and that set shrinks as we write matches.
    """
    # 1) Collect unique IDs to keep (first TAB-delimited field)
    keep_these = set()
    with open_text_maybe_gzip(chrom_file) as chromhandle:
        for line in chromhandle:
            line = line.strip()
            if not line:
                continue
            keep_these.add(line.split("\t", 1)[0])

    # 1) Open output (gz or plain)
    out_is_gz = output_fasta.endswith((".gz", ".gzip"))
    out = gzip.open(output_fasta, "wt") if out_is_gz else open(output_fasta, "w")
    write = out.write

    try:
        # 3) Stream over input FASTA using your parser; write matches and shrink the set
        for record in fasta.parse(input_fasta):
            pid = record.id
            if pid in keep_these:
                keep_these.discard(pid)     # let the set shrink to free RAM
                # Drop description to avoid later parsing issues; write directly (no format())
                write(f">{pid}\n")
                seq = record.seq            # parser already made this string
                # Write in chunks to avoid creating one giant formatted string
                for i in range(0, len(seq), 80):
                    write(seq[i:i+80] + "\n")
                if not keep_these:          # optional early exit if all found
                    break
    finally:
        out.close()

###### THESE ARE THE FUNCTIONS FOR ODP and ODP_SANDWICH
def genome_coords_to_plotstart_dict(path_to_genocoords_file, **kwargs):
    """
    Takes a genome coords file where:
      - col1: scaf name
      - col2: scaflen
      - col3: cumsum of the total plot size
      - col4: the plot starting position for that scaffold

    sca1 3822568 3822568  0
    sca2 2667796 6490364  3822568
    sca3 2526311 9016675  2667796
    sca4 2410750 11427425 2526311
    sca5 2150379 13577804 2410750
    sca6 1771964 15349768 2150379

    And returns a dict where col1 (scaf name) is key
     and col4 (plotting offset) is the value
    """
    offset_dict = {}
    with open(path_to_genocoords_file, "r") as f:
        for line in f:
            line = line.strip()
            if line:
                splitd = line.split()
                offset_dict[splitd[0]] = int(splitd[3])
    return offset_dict

def genome_coords_to_offset_dict(path_to_genocoords_file, **kwargs):
    """
    Takes a genome coords file where:
      - col1: scaf name
      - col2: scaflen
      - col3: cumsum of the total plot size
      - col4: the plot starting position for that scaffold

    sca1 3822568 3822568  0
    sca2 2667796 6490364  3822568
    sca3 2526311 9016675  2667796
    sca4 2410750 11427425 2526311
    sca5 2150379 13577804 2410750
    sca6 1771964 15349768 2150379

    And returns a dict where col1 (scaf name) is key
     and col3 (plotting offset) is the value
    """
    offset_dict = {}
    with open(path_to_genocoords_file, "r") as f:
        for line in f:
            line = line.strip()
            if line:
                splitd = line.split()
                offset_dict[splitd[0]] = int(splitd[2])
    return offset_dict

def generate_coord_structs_from_chrom_to_loc(chrom_file, **kwargs):
    """
    This parses a .chrom file and outputs five data structures that are easily
     used for mapping pandas dataframes.
    The output is a dict of dicts. Not the most intuitive format but easy for
     mapping to column values.
     { "prot_to_scaf":   prot_to_scaf,
       "prot_to_strand": prot_to_strand,
       "prot_to_start":  prot_to_start,
       "prot_to_stop":   prot_to_stop,
       "prot_to_middle": prot_to_middle }
    """
    prot_to_scaf   = {}
    prot_to_strand = {}
    prot_to_start  = {}
    prot_to_stop   = {}
    prot_to_middle = {}
    print("chrom_file", chrom_file)
    with open(chrom_file, "r") as f:
       for line in f:
           line = line.strip()
           if line:
               splitd = line.split()
               prot = splitd[0]
               # add things now
               prot_to_scaf[prot]   = splitd[1]
               prot_to_strand[prot] = splitd[2]
               start = int(splitd[3])
               prot_to_start[prot]  = start
               stop = int(splitd[4])
               prot_to_stop[prot]   = stop
               stop = int(splitd[4])
               prot_to_middle[prot] = int(start + (stop - start)/2)
    return { "prot_to_scaf":   prot_to_scaf,
             "prot_to_strand": prot_to_strand,
             "prot_to_start":  prot_to_start,
             "prot_to_stop":   prot_to_stop,
             "prot_to_middle": prot_to_middle }

def blast_plot_order_helper(coords, sample, xory, xprottoloc, yprottoloc, recip,
                            xorder, **kwargs):
    """
    This uses the reciprocal blast results to come up with the sort order
     for the y-axis scaffolds. Returns a list of the plot order.

    This code is all duplicated from the synteny plot function.
     Could be programmed in a better way to avoid redundancy, but this just fits
     the edge case where the y-axis has to be arranged based on the blast results.
    """
    # now make a lookup table of where the prots are.
    #  Use the x_offset and y_offset to recalculate where the plotting
    #  value is
    xcoords = generate_coord_structs_from_chrom_to_loc(xprottoloc)
    ycoords = generate_coord_structs_from_chrom_to_loc(yprottoloc)

    # now open the blast results and translate the pairs
    #  into plotting positions
    df = pd.read_csv(recip, header=None, sep = "\t")
    df.columns = ["xgene", "ygene", "pident", "length",
                  "mismatch", "gapopen", "qstart", "qend",
                  "sstart", "send", "evalue", "bitscore"]
    df = df[["xgene", "ygene", "bitscore", "evalue"]]

    #print(x_chrom)
    df["xpos"] = df["xgene"].map(xcoords["prot_to_middle"])
    df["ypos"] = df["ygene"].map(ycoords["prot_to_middle"])

    df["xscaf"] = df["xgene"].map(xcoords["prot_to_scaf"])
    df["yscaf"] = df["ygene"].map(ycoords["prot_to_scaf"])
    df = df.dropna()
    df = df.sort_values(by=['xpos'])
    df = df.dropna()

    grouped_df = df.groupby(["yscaf"])
    for key, item in grouped_df:
        max_item = grouped_df.get_group(key)['xscaf'].value_counts().idxmax()
        all_other_things = [x for x in grouped_df.get_group(key)['xscaf'].unique() if x != max_item]
        for thisthing in all_other_things:
            df = df.loc[~( (df["yscaf"] == key) & (df["xscaf"] == thisthing)), ]
    # now sort based on the xscafs and the xpos
    sorterIndex = dict(zip(xorder, range(len(xorder))))
    df.sort_values(['yscaf', 'ypos'],
        ascending = [True, True], inplace = True)
    df.reset_index(drop=True, inplace = True)
    df = df.drop_duplicates(subset=['yscaf'])
    df['x_Rank'] = df['xscaf'].map(sorterIndex)
    df.sort_values(['x_Rank', 'xpos'],
        ascending = [True, True], inplace = True)
    df = df.dropna()
    df.reset_index(drop=True, inplace = True)
    #print(list(df.yscaf))
    return(list(df.yscaf))

def parse_coords(coords_file, sample, xory,
                 xprottoloc=None, yprottoloc=None,
                 recip=None, xorder=None, **kwargs):
    """
    This parses the coordinates and returns a
      - coord-to-offset dict (I don't remember what this is for),
      - the size of each scaffold (a dictionary)
      - a list of locations to plot lines (These are the scaf/chrom divisions)
      - the max value for that axis
      - the tick labels
      - the tick positions
      - the yorder or xorder
    """
    config = kwargs["config"]
    offset = {}
    max_coord = 0
    lines_at = []
    df = pd.read_csv(coords_file, header = None, sep = " ")
    df.columns = ["scaf", "scaflen", "cumsum", "coordstart"]

    # if plotorder is in the config, then drop rows that aren't in it
    if "plotorder" in config["{}axisspecies".format(xory)][sample]:
        df = df[df["scaf"].isin(config["{}axisspecies".format(xory)][sample]["plotorder"])]
    # only drop if we haven't specified the plot order in the config
    #if drop_nas:
    #    df = df.dropna()
    df.reset_index(drop=True, inplace = True)
    df["cumsum"] = df["scaflen"].cumsum()
    df["cumsum"] = df["cumsum"] - df["scaflen"]
    print("df after cumulative sum and sorting")
    print(df)
    for i, row in df.iterrows():
        offset[row["scaf"]] = row["cumsum"]
        if i > 0:
            lines_at.append(row["cumsum"])
    max_coord = list(df["scaflen"].cumsum())[-1]

    #tick labels
    tick_labels = list(df["scaf"])
    tick_pos    = list(df["cumsum"] + (df["scaflen"]/2))

    scaf_to_len = {}
    for i, row in df.iterrows():
        scaf_to_len[row["scaf"]] = row["scaflen"]

    return (offset, scaf_to_len, lines_at, max_coord, tick_labels, tick_pos, list(df["scaf"]))



def determine_breaks(df, scaf_to_breaks_set, scaf_to_offset_dict,
                     sort_direction, auto_breaks, **kwargs):
    """
    determines the major breaks in Dx or Dy to use as partitions.

    The input parameters are:
      - df: the analysis df at the end of synteny plot.
      - sort_direction: either "x" or "y"

    The output of this method is a dataframe that is just the rows of the input
     that are the breakpoints in the input df.
    """
    # MAGIC NUMBERS
    # set window to change how many genes on either side are considered when
    #  looking for peaks. A value of 20 means 20 on either side, so 41 genes total
    window = 20
    smallwindow = 5
    # set small_window to resolve nearby peaks from different datasources

    sort_order = {"x": {"pos": "xmiddle",
                        "end": "xstop",
                        "chrom": "xscaf",
                        "D": "Dx"},
                  "y": {"pos": "ymiddle",
                        "end": "ystop",
                        "chrom": "yscaf",
                        "D": "Dy"}}

    # sort the dataframe based on which axis we're looking at
    df = df.sort_values(by=[sort_order[sort_direction]["pos"]])
    df = df.reset_index(drop=True)

    # first, figure out the manual break positions in terms of the protein coordinates
    manual_breaks_indices = set()
    for thisscaf in scaf_to_breaks_set:
        for thisposition in scaf_to_breaks_set[thisscaf]:
            offset_position = thisposition + scaf_to_offset_dict[thisscaf]
            subdf = df.loc[df[sort_order[sort_direction]["end"]] <= offset_position, ]
            subdf = subdf.sort_values(by=[sort_order[sort_direction]["end"]])
            tempdf = df.loc[df[sort_order[sort_direction]["chrom"]] == thisscaf, ["xgene", "ygene", "xstart"]]
            #print(thisscaf, thisposition, offset_position)
            #print(tempdf)
            #print(subdf.loc[subdf.index[-1]])
            #print()
            #sys.exit()
            manual_breaks_indices.add(subdf.index[-1])
    manual_breaks_indices = list(manual_breaks_indices)
    #if not auto_breaks:
    #    print("manual_breaks")
    #    print(manual_breaks_indices)

    all_ranges = set()
    if auto_breaks:
        unique_chroms = []
        #chrom_breakpoints = []
        for index, row in df.iterrows():
            thischrom = row[sort_order[sort_direction]["chrom"]]
            thispos  = row[sort_order[sort_direction]["pos"]]
            if thischrom not in unique_chroms:
                unique_chroms.append(thischrom)
                #chrom_breakpoints.append(thispos)
        ## this line is solely for plotting. Not useful
        #chrom_breakpoints = chrom_breakpoints[1::]

        # There are three different analysis types that we will use to figure
        #  out the seps.
        # - deltMA is the derivative of the smoothed data
        # - deltD is the derivative of the raw data
        # - Dx2 is the raw D data.
        #
        # All of the data are selected based on the max value of what was above
        #  the median for that chromosome.
        for thiscol in ["MA", "deltMA","deltD", "Dx2"]:
            for thischrom in unique_chroms:
                # use .copy() to make sure we're not modifying the original df
                subdf = df.loc[df[sort_order[sort_direction]["chrom"]] == thischrom, ].copy()
                # Dx2 is just the raw data that is above the median
                subdf["Dx2"] = subdf[sort_order[sort_direction]["D"]]
                subdf['Dx2'] = np.where((subdf[sort_order[sort_direction]["D"]] < subdf[sort_order[sort_direction]["D"]].median()),np.nan,subdf["Dx2"])
                # MA is the moving average of the raw data
                subdf["MA"] = subdf["Dx2"].rolling(window=3, center=True).mean()
                subdf["MA2"] = subdf["Dx2"].rolling(window=19, center=True).mean()
                # deltMA is the derivative of the moving average
                subdf["deltMA"] = subdf["MA"].diff() / subdf["MA"].index.to_series().diff()
                subdf['deltMA'] = np.where((subdf["MA"] < subdf["MA"].median()),np.nan,subdf["deltMA"])
                # deltD is the derivative of the raw data
                subdf["deltD"] = subdf["Dx2"].diff() / subdf["Dx2"].index.to_series().diff()
                subdf['deltD'] = np.where((subdf.Dx2 < subdf["Dx2"].median()),np.nan,subdf.deltD)

                # get the groups of consecutive values in each category
                idxmaxes = set()
                ind = list(subdf[~subdf[thiscol].isnull()].index)
                ranges =[]
                for k,g in groupby(enumerate(ind),lambda x:x[0]-x[1]):
                    group = (map(itemgetter(1),g))
                    group = list(map(int,group))
                    ranges.append((group[0],group[-1]))

                # now get the peak from each contiguous range of values
                if len(ranges) > 0:
                    for this_range in ranges:
                        if this_range[0] != this_range[-1]:
                            #this_range = [x for x in range(this_range[0], this_range[1]+1)]
                            this_range = list(this_range)
                            which_d_col = sort_order[sort_direction]["D"]
                            temp = subdf.loc[this_range[0]:this_range[-1]][which_d_col].idxmax()
                            idxmaxes.add(temp)

                # picks the best in large windows of genes.
                #  See the description for the 'window' variable above
                keep_idx_maxes = set()
                ignore_set = set()
                done = False
                consider_ranges = set()
                while not done:
                    consider_ranges = set()
                    # get peaks within the window if they're not in the ignore set
                    for this_idx in idxmaxes:
                        thistup = tuple([x for x in idxmaxes
                             if ((x > this_idx - window)
                                 and (x < this_idx + window)
                                 and (x not in ignore_set))])
                        if len(thistup) > 0:
                            consider_ranges.add(thistup)

                    # now for each set of peaks, get the best in each window
                    consider_ranges = sorted(list(consider_ranges), key=len, reverse=True)
                    if len(consider_ranges) > 0: # skip the empty ranges
                        thisrange = list(consider_ranges[0])
                        if len(thisrange) == 1:
                            done = True
                        else:
                            submax = df.loc[thisrange, ][sort_order[sort_direction]["D"]].idxmax()
                            for thisid in thisrange:
                                if thisid != submax:
                                    ignore_set.add(thisid)
                    else: # if it is empty, leave
                        done = True
                # We found the biggest peaks in the windows, add them to all_ranges
                for entry in consider_ranges:
                    all_ranges.add(entry)

    # flatten the results of what we got from the last analysis
    idxmaxes = flatten(all_ranges)
    idxmaxes = flatten([idxmaxes, manual_breaks_indices])

    # From the dataset of all peaks, find the best in small windows.
    #  See the variable `smallwindow` above
    # The reason we have this block is that the same peak, or something near it,
    #  could have been added multiple times, at slightly different indices.
    #  This collapses the similar indices to get the best.
    ignore_set = set()
    done = False
    consider_ranges = set()
    while not done:
        consider_ranges = set()
        for this_idx in idxmaxes:
            # get windows of ranges if they're not in the ignore set
            thistup = tuple([x for x in idxmaxes
                 if ((x > this_idx - smallwindow)
                     and (x < this_idx + smallwindow)
                     and (x not in ignore_set))])
            if len(thistup) > 0:
                consider_ranges.add(thistup)

        consider_ranges = sorted(list(consider_ranges), key=len, reverse=True)
        if len(consider_ranges) > 0:
            thisrange = list(consider_ranges[0])
            if len(thisrange) == 1:
                done = True
            else:
                submax = df.loc[thisrange, ][sort_order[sort_direction]["D"]].idxmax()
                for thisid in thisrange:
                    if thisid != submax:
                        ignore_set.add(thisid)
        else:
            # there's nothing here
            done = True

    # vert_lines is the list of indices from the df that have the peaks that
    #  we want to keep.
    vert_lines = flatten(consider_ranges)
    # return a dataframe of the intersections we want
    return df.loc[vert_lines].copy()

