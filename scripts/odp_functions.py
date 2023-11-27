"""
These are functions that are shared by odp, odp_trio, and other scripts
"""
# this is all needed to load our custom fasta parser
import os
import sys
snakefile_path = os.path.dirname(os.path.abspath(__file__)) 
dependencies_path = os.path.join(snakefile_path, "../dependencies/fasta-parser")
sys.path.insert(1, dependencies_path)
import fasta

# ODP-specific imports
import odp_plotting_functions as odp_plot

# other standard python libraries
from itertools import combinations
from itertools import product

# non-standard dependencies
import pandas as pd

# TODO implement this function soon after release
#def gff_or_chrom(filepath):
#    """
#    This function takes a filepath and determines if it is a gff or a chrom file.
#    Returns "gff" or "chrom", or throws an error if the software can't decide.
#    """
#    # First check what the file extension is

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
    odp_path = os.path.abspath(os.path.join(snakefile_path, ".."))
    cwd      = os.getcwd()

    # test if we are in the odp directory
    if odp_path in cwd:
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
    fdf = (f_raw.groupby("qseqid")
             .apply(lambda group: group.loc[group["evalue"] == group['evalue'].min()])
             .reset_index(drop=True)
          )

    r_raw = pd.read_csv(y_to_x_blastp_results, sep="\t")
    r_raw.columns = ["qseqid", "sseqid", "pident", "length",
                   "mismatch", "gapopen", "qstart", "qend",
                   "sstart", "send", "evalue", "bitscore"]
    rdf = (r_raw.groupby("qseqid")
             .apply(lambda group: group.loc[group["evalue"] == group['evalue'].min()])
             .reset_index(drop=True)
          )
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
    fdf = f_raw.copy()
    fdf = (fdf.groupby("qseqid")
             .apply(lambda group: group.loc[group["evalue"] == group['evalue'].min()])
             .reset_index(drop=True)
          )
    rdf = r_raw.copy()
    rdf = (rdf.groupby("qseqid")
             .apply(lambda group: group.loc[group["evalue"] == group['evalue'].min()])
             .reset_index(drop=True)
          )

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

def reciprocal_best_hits_jackhmmer(
        x_to_y_blastp_results, y_to_x_blastp_results, outfile):
    """
    This function finds reciprocal best jackhmmer hits between two samples.
    The input is a blastp results file where x was jackhmmer'd against y,
      and a blastp results file where y was jackhmmer'd against x.

    The output format is just the rows of the blastp results from the x_to_y file.
    Saves it as a df to outfile.

    This algorithm is permissive in that it finds the best hits between the two
      species even if the e-values for the "best hit" are equivalent. This fixes
      one of the problems with blastp results. The results are still reciprocal
      best, though.
    """
    jackhmmercol = ["target_name", "accession",  "query_name",
                    "accession",
                    "evalue",  "score",          "bias",
                    "dom_evalue2", "dom_score2", "bias2",
                    "exp", "reg", "clu",
                    "ov", "env", "dom", "rep", "inc",
                    "description_of_target"]
    f_raw = pd.read_csv(x_to_y_blastp_results,
                        sep = "\s+", comment = "#",
                        usecols=range(len(jackhmmercol)))
    f_raw.columns = jackhmmercol
    fdf = f_raw.sort_values(["query_name", "score", "evalue" ], ascending=[True, False, True]).drop_duplicates(subset="query_name")
    #fdf = f_raw.sort_values(["query_name", "score", "evalue" ], ascending=[True, False, True]).groupby("query_name").head(2)


    r_raw = pd.read_csv(y_to_x_blastp_results,
                        sep="\s+", comment = "#",
                        usecols=range(len(jackhmmercol)))
    r_raw.columns = jackhmmercol
    rdf = r_raw.sort_values(["query_name", "score", "evalue"], ascending=[True, False, True]).drop_duplicates(subset="query_name")
    #rdf = r_raw.sort_values(["query_name", "score", "evalue"], ascending=[True, False, True]).groupby("query_name").head(2)

    rdf.columns = ["query_name", "accession",  "target_name",
                    "accession",
                    "evalue",  "score",          "bias",
                    "dom_evalue2", "dom_score2", "bias2",
                    "exp", "reg", "clu",
                    "ov", "env", "dom", "rep", "inc",
                    "description_of_target"]

    rdf = rdf[["target_name","query_name"]]

    #These are the singleton RBH
    new_df = pd.merge(fdf, rdf,  how='inner',
                      left_on  = ['query_name','target_name'],
                      right_on = ['query_name','target_name'])
    #these rows are a little pedantic and we don't really need to do them
    new_df = new_df.sort_values(["query_name","score"],
                                ascending=[True, False]).drop_duplicates(
                                    subset="query_name")
    new_df = new_df.sort_values(["target_name","score"],
                                ascending=[True, False]).drop_duplicates(
                                    subset="query_name")

    new_df.columns = ["sseqid", "accession",  "qseqid",
                      "accession",
                      "evalue",  "bitscore",          "bias",
                      "dom_evalue2", "dom_score2", "bias2",
                      "exp", "reg", "clu",
                      "ov", "env", "dom", "rep", "inc",
                      "description_of_target"]
    new_df["pident"]   = 0
    new_df["length"]   = 0
    new_df["mismatch"] = 0
    new_df["gapopen"]  = 0
    new_df["qstart"]   = 0
    new_df["qend"]     = 0
    new_df["sstart"]   = 0
    new_df["send"]     = 0
    new_df = new_df[["qseqid", "sseqid", "pident", "length",
                   "mismatch", "gapopen", "qstart", "qend",
                   "sstart", "send", "evalue", "bitscore"]]
    print(new_df)
    new_df.to_csv(outfile, sep="\t", index = False, header = False)

def check_file_exists(filepath) -> bool:
    """
    checks if a file exists.
    If not, raises an error
    """
    if not os.path.isfile(filepath):
        outmessage =  "*********************************************************************\n"
        outmessage += "* ERROR:\n"
        outmessage += "*  This file does not exist:" + filepath + "\n"
        outmessage =  "*********************************************************************\n"
        raise IOError(outmessage)
    else:
        return True

def chrom_file_is_legal(chrompath):
    """
    Checks if a chrom file is legal.
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
    # go through the file line by line and inspect each element
    with open(chrompath, "r") as f:
        for line in f:
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


def check_species_input_legality(fastapath, peppath, chrompath) -> bool:
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
        outmessage += "*  There are " + str(len(duplicate)) + " duplicate sequence headers.\n"
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
    protein_headers     = set()
    duplicate_headers   = set()
    protein_sequences   = set()
    duplicate_sequences = set()
    for record in fasta.parse(peppath):
        if record.id not in protein_headers:
            protein_headers.add(record.id)
        else:
            duplicate_headers.add(record.id)
        if str(record.seq) not in protein_sequences:
            protein_sequences.add(str(record.seq))
        else:
            duplicate_sequences.add(record.id)

    if len(duplicate_headers) > 0:
        dupstring = "".join(["*    - " + str(x) + "\n" for x in sorted(duplicate_headers)[:3]])
        # raise an error because each ID should only occur once
        outmessage =  "*********************************************************************\n"
        outmessage += "* ERROR:\n"
        outmessage += "*  There is a protein fasta with duplicate sequence headers.\n"
        outmessage += "*  Each sequence in the protein fasta must have a unique ID.\n"
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

    # 3. Check that each sequence ID exists only once
    if len(duplicate_sequences) > 0:
        dupstring = "".join(["*    - " + str(x) + "\n" for x in sorted(duplicate_sequences)[:3]])
        # raise an error because each ID should only occur once
        outmessage =  "*********************************************************************\n"
        outmessage += "* ERROR:\n"
        outmessage += "*  Some protein sequences in your file are identical.\n"
        outmessage += "*  Each protein sequence must be unique.\n"
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
        outmessage += "*********************************************************************\n"
        raise IOError(outmessage)

    # PARSE AND CHECK THE CHROM FILE
    # 1. check that the file exists
    check_file_exists(chrompath)
    proteins_not_in_pep    = set()
    scaffolds_not_in_fasta = set()
    for line in open(chrompath, 'r'):
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
        # raise an error because the proteins should have been seen already
        dupstring = "".join(["*    - " + str(x) + "\n" for x in sorted(proteins_not_in_pep)[:3]])
        outmessage =  "*********************************************************************\n"
        outmessage += "* ERROR:\n"
        outmessage += "*  Some proteins in the .chrom file were not seen in the protein\n"
        outmessage += "*   .fasta file.\n"
        outmessage += "*\n"
        outmessage += "*  The chrom file with the problem is: " + chrompath + "\n"
        outmessage += "*  There are " + str(len(proteins_not_in_pep)) + " proteins in the .chrom not seen in the protein .fasta\n"
        outmessage += "*  Here are the first 1 to 3:\n"
        outmessage += dupstring
        outmessage += "*\n"
        outmessage += "*  The reason this is problematic is that we need to access every\n"
        outmessage += "*   protein specified in the .chrom file, but it is unavailable.\n"
        outmessage += "*\n"
        outmessage += "*  Please investigate whether there are too many entries in the .chrom\n"
        outmessage += "*   file, or if something is missing from the protein .fasta file.\n"
        outmessage += "*   Then, fix your files and re-run this pipeline.\n"  
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
    legal = ["proteins", "chrom", "genome", "genus",
             "minscafsize", "manual_breaks", "chrom_to_color",
             "plotorder", "species", "prot_to_group"]
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


def generate_coord_structs_from_chrom_to_loc(chrom_file):
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

def filter_fasta_chrom(chrom_file, input_fasta, output_fasta):
    """
    takes a chrom file, only keeps proteins in input_fasta from chrom file,
     saves those prots to output_fasta
    """
    keep_these = set()
    printed_already = set()
    with open(chrom_file, "r") as f:
        for line in f:
            line = line.strip()
            if line:
                splitd = line.split()
                keep_these.add(splitd[0])
    outhandle = open(output_fasta, "w")
    for record in fasta.parse(input_fasta):
        if record.id in keep_these and record.id not in printed_already:
            # The record object has the properties
            #  - Record.id
            #  - Record.seq
            #  - Record.desc
            # get rid of the description to avoid parsing errors later
            record.desc=""
            print(record, file = outhandle)
            printed_already.add(record.id)
    outhandle.close()

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
    df["yscaf"] = df["ygene"].map(xcoords["prot_to_scaf"])
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
    # now figure out if we need to sort or not
    drop_nas = True

    plotorder = None

    #print("df after xory")
    # now we have determined if we need to sort
    if plotorder != None: # if plotorder has something in it.
        #print(" - using custom plot order: ", plotorder)
        sortdict = {key: val for key, val in zip(plotorder, range(len(plotorder)))}
        df['rank'] = df['scaf'].map(sortdict)
        df.sort_values(by = 'rank' ,inplace=True)
    #print("df after plotorder")
    #print(df)

    # now, if we made plotorder from config then drop rows
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


def calc_D_for_y_and_x(df, x_offset, y_offset, x_scaf_to_len, y_scaf_to_len, **kwargs):
    """
    This calculates D for both the x and y axes.
    Defined in the 2020 vertebrate synteny paper.
    """
    df = df.dropna()
    # some variable names in this for loop are "x" but it doesn't matter.
    #  everything important is variable between x and y
    for thisdir in ["x", "y"]:
        df = df.sort_values(by=["{}middle".format(thisdir)])
        df.reset_index(drop=True, inplace = True)

        unique_x = df["{}scaf".format(thisdir)].unique()
        thisdir_dfs = []
        # this just calculates Dx
        for thisx in unique_x:
            xdf = df.loc[df["{}scaf".format(thisdir)] == thisx, ].copy()
            xdf = xdf.reset_index(drop=True)
            oppositexy = "d"
            this_offset = {}
            this_scaf_to_len = {}
            if thisdir == "x":
                oppositexy = "y"
                this_offset = x_offset
                this_scaf_to_len = x_scaf_to_len
            elif thisdir == "y":
                oppositexy = "x"
                this_offset = y_offset
                this_scaf_to_len = y_scaf_to_len
            df2 = pd.get_dummies(xdf["{}scaf".format(oppositexy)])
            df2_xiL = df2.apply(lambda x: x.rolling(20).mean(), axis = 0)
            df2_xiR = df2.apply(lambda x: x.iloc[::-1].rolling(20).mean(), axis = 0).iloc[::-1]
            df2_xiR = df2_xiR.set_index(df2_xiR.index - 1)
            df2_xiR = df2_xiR.iloc[1:]
            subtractdf = df2_xiR.fillna(0) - df2_xiL.fillna(0)
            D = subtractdf.apply(lambda x: np.sqrt(np.square(x).sum()), axis = 1)
            xdf["D{}".format(thisdir)] = D
            xdf["D{}_barleft".format(thisdir)] = 0
            xdf["D{}_barmiddle".format(thisdir)] = 0
            xdf["D{}_barright".format(thisdir)] = 0
            xdf["D{}_barwidth".format(thisdir)] = 0
            for i, row in xdf.iterrows():
                barleft   = -1
                barright  = -1
                barmiddle = -1
                barwidth  = -1
                if len(xdf) > 1:
                    if i == 0:
                        thisend   = row["{}stop".format(thisdir)]
                        nextstart = xdf.loc[i+1, "{}start".format(thisdir)]
                        barleft   = this_offset[thisx]
                        barright  = thisend + ((nextstart-thisend)/2)
                    elif i == (len(xdf) - 1):
                        prevend   = xdf.loc[i-1, "{}stop".format(thisdir)]
                        thisstart = row["{}start".format(thisdir)]
                        barleft   = prevend + ((thisstart-prevend)/2)
                        barright  = this_scaf_to_len[thisx]
                    else:
                        prevend   = xdf.loc[i-1, "{}stop".format(thisdir)]
                        thisstart = row["{}start".format(thisdir)]
                        thisend   = row["{}stop".format(thisdir)]
                        nextstart = xdf.loc[i+1, "{}start".format(thisdir)]
                        barleft   = prevend + ((thisstart-prevend)/2)
                        barright  = thisend + ((nextstart-thisend)/2)
                xdf.loc[i, "D{}_barleft".format(thisdir)]   = barleft
                xdf.loc[i, "D{}_barright".format(thisdir)]  = barright
                xdf.loc[i, "D{}_barmiddle".format(thisdir)] = barleft + ((barright - barleft)/2)
                xdf.loc[i, "D{}_barwidth".format(thisdir)]  = barright - barleft + 1
            thisdir_dfs.append(xdf)
        df = pd.concat(thisdir_dfs)
    df = df.sort_values(by=["xmiddle"])
    df.reset_index(drop=True, inplace = True)
    return df

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
                subdf['Dx2'] = np.where((subdf[sort_order[sort_direction]["D"]] < subdf[sort_order[sort_direction]["D"]].median()),np.NaN,subdf["Dx2"])
                # MA is the moving average of the raw data
                subdf["MA"] = subdf["Dx2"].rolling(window=3, center=True).mean()
                subdf["MA2"] = subdf["Dx2"].rolling(window=19, center=True).mean()
                # deltMA is the derivative of the moving average
                subdf["deltMA"] = subdf["MA"].diff() / subdf["MA"].index.to_series().diff()
                subdf['deltMA'] = np.where((subdf["MA"] < subdf["MA"].median()),np.NaN,subdf["deltMA"])
                # deltD is the derivative of the raw data
                subdf["deltD"] = subdf["Dx2"].diff() / subdf["Dx2"].index.to_series().diff()
                subdf['deltD'] = np.where((subdf.Dx2 < subdf["Dx2"].median()),np.NaN,subdf.deltD)

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

def gen_plotting_df(ycoords_file, xcoords_file,
                    xprottoloc, yprottoloc,
                    xsample, ysample,
                    recip, outtable, plotorder_file, **kwargs):
    """
    Generates a dataframe that will be used by the other parts of the program
     for plotting.
    Saves it to a file
    """
    import pandas as pd
    import numpy as np

    # first make a lookup table of how to calculate the
    #  x and y coords_file. This lookup is just the amount of
    # bp to add to the value when plotting. We pass the xchrom,
    #  xprot_to_scaf in case we need to sort everything based on order of
    #  occurrence on the scaffolds
    x_offset, x_scaf_to_len, vertical_lines_at, xmax, xticklabel, xtickpos, xorder = parse_coords(
        xcoords_file, xsample, "x", **kwargs)
    print("found {} x chromosomes".format(len(x_offset)))

    y_offset, y_scaf_to_len, horizontal_lines_at, ymax, yticklabel, ytickpos, yorder = parse_coords(
        ycoords_file, ysample, "y",
        xprottoloc, yprottoloc, recip, xticklabel)
    print("found {} y chromosomes".format(len(y_offset)))

    # now save the plot order to a file
    with open(plotorder_file, "w") as f:
        print("xplotorder:", file=f)
        for entry in xorder:
            print("  - {}".format(entry), file=f)
        print("yplotorder:", file=f)
        for entry in yorder:
            print("  - {}".format(entry), file=f)

    # now make a lookup table of where the prots are.
    #  Use the x_offset and y_offset to recalculate where the plotting
    #  value is
    xstruct = generate_coord_structs_from_chrom_to_loc(xprottoloc)
    xstruct["prot_plot_start"] = {}
    xstruct["prot_plot_middle"] = {}
    xstruct["prot_plot_stop"] = {}
    ystruct = generate_coord_structs_from_chrom_to_loc(yprottoloc)
    ystruct["prot_plot_start"] = {}
    ystruct["prot_plot_middle"] = {}
    ystruct["prot_plot_stop"] = {}
    #print("xstruct")
    #print(xstruct)
    # get rid of proteins that we don't need along the x axis
    #  also set the plotting position based on the offset
    for prot in list(xstruct["prot_to_middle"].keys()):
        scaf = xstruct["prot_to_scaf"][prot]
        if scaf not in x_offset:
            for thisdict in xstruct:
                xstruct[thisdict].pop(prot, None)
        else:
            xstruct["prot_plot_start"][prot]  = xstruct["prot_to_start"][prot] + x_offset[scaf]
            xstruct["prot_plot_middle"][prot] = xstruct["prot_to_middle"][prot] + x_offset[scaf]
            xstruct["prot_plot_stop"][prot]   = xstruct["prot_to_stop"][prot] + x_offset[scaf]
    # get rid of proteins that we don't need along the y axis
    #  also set the plotting position based on the offset
    for prot in list(ystruct["prot_to_middle"].keys()):
        scaf = ystruct["prot_to_scaf"][prot]
        if scaf not in y_offset:
            for thisdict in xstruct:
                ystruct[thisdict].pop(prot, None)
        else:
            ystruct["prot_plot_start"][prot]  = ystruct["prot_to_start"][prot] + y_offset[scaf]
            ystruct["prot_plot_middle"][prot] = ystruct["prot_to_middle"][prot] + y_offset[scaf]
            ystruct["prot_plot_stop"][prot]   = ystruct["prot_to_stop"][prot] + y_offset[scaf]

    # now open the blast results and translate the pairs
    #  into plotting positions
    df = pd.read_csv(recip, header=None, sep = "\t")
    df.columns = ["xgene", "ygene", "pident", "length",
                  "mismatch", "gapopen", "qstart", "qend",
                  "sstart", "send", "evalue", "bitscore"]
    df = df[["xgene", "ygene", "bitscore", "evalue"]]
    df["xgene"] = df["xgene"].astype(str)
    df["ygene"] = df["ygene"].astype(str)
    #print(x_chrom)
    df["xstart"]  = df["xgene"].map(xstruct["prot_plot_start"])
    df["xmiddle"] = df["xgene"].map(xstruct["prot_plot_middle"])
    df["xstop"]   = df["xgene"].map(xstruct["prot_plot_stop"])

    df["ystart"]  = df["ygene"].map(ystruct["prot_plot_start"])
    df["ymiddle"] = df["ygene"].map(ystruct["prot_plot_middle"])
    df["ystop"]   = df["ygene"].map(ystruct["prot_plot_stop"])
    print(df)

    # I don't remember why the x and y is switched here.
    df["yscaf"] = df["ygene"].map(ystruct["prot_to_scaf"])
    df["xscaf"] = df["xgene"].map(xstruct["prot_to_scaf"])
    #df = df.dropna() # this messed up plotting order if there were no hits on that scaf.
    df = df.sort_values(by=['xmiddle'])
    #df = df.loc[df["evalue"] <= float("1E-20"), ]
    df.reset_index(drop=True, inplace = True)

    # Now calculate D and Dw for both X and Y axes
    df = calc_D_for_y_and_x(df, x_offset, y_offset, x_scaf_to_len, y_scaf_to_len)
    if outtable:
        df.to_csv(outtable, sep="\t")

def synteny_plot(plotting_df,    xcoords_file,  ycoords_file,
                 xsample,        ysample,
                 xbreaks_file,   ybreaks_file,
                 synplot,        prot_to_color,  dropmissing,
                 plot_x_lines = False,
                 plot_y_lines = False,
                 xprottoloc = False,
                 yprottoloc = False, **kwargs):
    """
    If the user provided a plot order, then we should not skip any scaffolds.

    This is the main plotting script for the synteny plot
    """
    config = kwargs["config"]
    import pandas as pd
    import matplotlib
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    #import matplotlib.patches as mplpatches
    from matplotlib.ticker import StrMethodFormatter, NullFormatter
    import numpy as np

    # CALL THIS TO GET THE VISUAL STYLE WE NEED
    odp_plot.format_matplotlib()

    # Preserve the vertical order of embedded images:
    matplotlib.rcParams['image.composite_image'] = False
    # text as font in pdf
    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42

    # open the xbreaks and the ybreaks dataframes
    # These are where there are low-opacity,
    #  dotted lines to show breaks in synteny
    # example entries are like this:
    xbreaks_df = pd.read_csv(xbreaks_file, delimiter="\t", index_col=0)
    xbreaks_df["xgene"] = xbreaks_df["xgene"].astype(str)
    xbreaks_df["ygene"] = xbreaks_df["ygene"].astype(str)
    ybreaks_df = pd.read_csv(ybreaks_file, delimiter="\t", index_col=0)
    ybreaks_df["xgene"] = ybreaks_df["xgene"].astype(str)
    ybreaks_df["ygene"] = ybreaks_df["ygene"].astype(str)


    # first make a lookup table of how to calculate the
    #  x and y coords. This lookup is just the amount of
    # bp to add to the value when plotting. We pass the xchrom,
    #  xprot_to_scaf in case we need to sort everything based on order of
    #  occurrence on the scaffolds
    x_offset, x_scaf_to_len, vertical_lines_at, xmax, xticklabel, xtickpos, xorder = parse_coords(
        xcoords_file, xsample, "x", **kwargs)
    print("found {} x chromosomes".format(len(x_offset)))

    y_offset, y_scaf_to_len, horizontal_lines_at, ymax, yticklabel, ytickpos, yorder = parse_coords(
        ycoords_file, ysample, "y")
    print("found {} y chromosomes".format(len(y_offset)))

    # first make a lookup table
    df = pd.read_csv(plotting_df, delimiter="\t", index_col=0)
    df["xgene"] = df["xgene"].astype(str)
    df["ygene"] = df["ygene"].astype(str)

    xgene = df["xgene"]
    x = df["xmiddle"]
    y = df["ymiddle"]
    bitscore = df["bitscore"]

    print("found {} points to plot".format(len(df["xmiddle"])))
    print("max bitscore: ", max(bitscore))
    bitscore_adjusted = [(x/max(bitscore))*60 for x in bitscore]
    colors = [(0, 0, 1.0, min( 1.0, (x/max(bitscore))*(max(bitscore)/np.mean(bitscore)))) for x in bitscore]
    drops = set()
    if prot_to_color:
        for i in range(len(xgene)):
            alpha = colors[i][3]
            try:
                newcolor = list(matplotlib.colors.to_rgba(prot_to_color[xgene[i]]))
            except:
                print("couldn't find a color for: {}".format(xgene[i]))
                newcolor = [0,0,0,1]
                drops.add(i)
            newcolor[3] = alpha
            colors[i] = newcolor

    if dropmissing:
        vararray = [xgene, x, y, bitscore, bitscore_adjusted, colors]
        for j in range(len(vararray)):
            vararray[j] = [vararray[j][i] for i in range(len(vararray[j])) if i not in drops]
        xgene = vararray[0]
        x     = vararray[1]
        y     = vararray[2]
        bitscore = vararray[3]
        bitscore_adjusted = vararray[4]
        colors = vararray[5]

    # now make a scatter plot
    figWidth = 8
    figHeight = 8
    plt.figure(figsize=(figWidth,figHeight))
    #set the panel dimensions
    panelWidth = 4
    panelHeight = 4
    dpanel_width = 0.25
    #find the margins to center the panel in figure
    leftMargin = (figWidth - panelWidth)/2
    bottomMargin = ((figHeight - panelHeight)/2)
    panel1 = plt.axes([leftMargin/figWidth, #left
                         bottomMargin/figHeight,    #bottom
                         panelWidth/figWidth,   #width
                         panelHeight/figHeight])     #height
    panelxd = plt.axes([leftMargin/figWidth, #left
                         (bottomMargin+panelHeight+0.1)/figHeight,    #bottom
                         panelWidth/figWidth,   #width
                         dpanel_width/figHeight])     #height
    panelyd = plt.axes([(leftMargin+panelWidth + 0.1)/figWidth, #left
                         bottomMargin/figHeight,    #bottom
                         dpanel_width/figWidth,   #width
                         panelHeight/figHeight])     #height
    panel1.tick_params(axis='both',which='both',
                        bottom=False, labelbottom=True,
                        left=False, labelleft=True,
                        right=False, labelright=False,
                        top=False, labeltop=False)
    panelxd.tick_params(axis='both',which='both',
                        bottom=False, labelbottom=False,
                        left=False, labelleft=False,
                        right=False, labelright=False,
                        top=False, labeltop=False)
    panelyd.tick_params(axis='both',which='both',
                        bottom=False, labelbottom=False,
                        left=False, labelleft=False,
                        right=False, labelright=False,
                        top=False, labeltop=False)
    # set the panel linewidth thinner
    for this_panel in [panel1, panelxd, panelyd]:
        for axis in ['top','bottom','left','right']:
            this_panel.spines[axis].set_linewidth(0.5)
    # turn off the axis spines
    for this_panel in [panelxd, panelyd]:
        this_panel.spines['top'].set_visible(False)
        this_panel.spines['right'].set_visible(False)

    panel1.scatter(x, y, color = colors,
                   ec = None, s=6, linewidths = 0)
    print("xmax is")
    print(xmax)
    # set mins and max
    panel1.set_xlim([0, xmax])
    panel1.set_ylim([0, ymax])

    # set x ticks
    newarr = []
    newarrlabels=[]
    if not plot_x_lines:
        #there are inevitably going to be many scaffolds. We need to subset
        # get a list of evenly spaced indices
        numElems = min(20, len(xtickpos)) # this could break if there are fewer elements
        arr = xtickpos
        idx = np.round(np.linspace(0, len(arr) - 1, numElems)).astype(int)
        newarr       = [arr[i] for i in idx]
        newarrNumScaf     = [i for i in idx]
        # turn on y-axis ticks on the left - plot scaffolds
        panel1.tick_params(bottom=True)
        panel1.set_xticks(newarr)
        panel1.set_xticklabels(newarrNumScaf, fontsize=8, rotation=90)
        panel1.set_xlabel(xsample + " number of scaffolds")
    else:
        newarr = [0] + list(x_offset.values())
        panel1.set_xticks(xtickpos)
        panel1.set_xticklabels(xticklabel, fontsize=8, rotation = 90)
        panel1.set_xlabel(xsample + " scaffolds")
    # turn on x-axis ticks on the Dx plot
    newarrlabels = [round(x/1000000, 1) for x in newarr]
    panelxd.tick_params(top=True, labeltop=True)
    panelxd.set_xticks(newarr)
    panelxd.set_xticklabels(newarrlabels, fontsize=8, rotation=90)
    panelxd.xaxis.set_label_position("top")
    panelxd.set_xlabel("Mb")

    # set y ticks
    newarr=[]
    newarrlabels=[]
    if not plot_y_lines:
        #there are inevitably going to be many scaffolds. We need to subset
        # get a list of evenly spaced indices
        numElems = min(20, len(ytickpos))
        arr = ytickpos
        idx = np.round(np.linspace(0, len(arr) - 1, numElems)).astype(int)
        newarr       = [arr[i] for i in idx]
        newarrlabels = [round(arr[i]/1000000, 1) for i in idx]
        newarrNumScaf     = [i for i in idx]
        # turn on y-axis ticks on the left - plot scaffolds
        panel1.tick_params(left=True)
        panel1.set_yticks(newarr)
        panel1.set_yticklabels(newarrNumScaf, fontsize=8)
        panel1.set_ylabel(ysample + " number of scaffolds")
    else:
        newarr = [0] + horizontal_lines_at + [ymax]
        panel1.set_yticks(ytickpos)
        panel1.set_yticklabels(yticklabel, fontsize=8)
        panel1.set_ylabel(ysample + " scaffolds")
    # turn on y-axis ticks on the Dy plot
    newarrlabels = [round(x/1000000, 1) for x in newarr]
    panelyd.tick_params(right=True, labelright=True)
    panelyd.set_yticks(newarr)
    panelyd.set_yticklabels(newarrlabels, fontsize=8)
    panelyd.yaxis.set_label_position("right")
    panelyd.set_ylabel("Mb")

    # set the x and y labels on Dy and Dx
    panelxd.bar(x = df["Dx_barmiddle"], height=df["Dx"], width = df["Dx_barwidth"],
                align = "center", lw=0, color="blue", zorder = 2)
    panelxd.set_xlim([0,xmax])
    panelxd.set_ylabel('Dx', fontsize=10)

    panelyd.barh(y = df["Dy_barmiddle"], width=df["Dy"], height = df["Dy_barwidth"],
                 align = "center", lw=0, color="blue", zorder = 2)
    panelyd.set_ylim([0,ymax])
    panelyd.set_xlabel('Dy', fontsize=10)

    for this_axis in [panel1, panelxd, panelyd]:
        this_axis.xaxis.get_offset_text().set_visible(False)
        this_axis.yaxis.get_offset_text().set_visible(False)

    #plot vertical lines
    if plot_x_lines:
        for value in vertical_lines_at:
            panel1.axvline(x=value, color="black", lw=0.5)
    #plot horizontal lines
    if plot_y_lines:
        for value in horizontal_lines_at:
            panel1.axhline(y=value, color="black", lw=0.5)

    # plot vertical BOS
    for value in xbreaks_df["xmiddle"]:
        panel1.axvline(x=value, color=[0,0,0,0.25], lw=0.5, linestyle="dotted")
    # plot horizontal BOS
    for value in ybreaks_df["ymiddle"]:
        panel1.axhline(y=value, color=[0,0,0,0.25], lw=0.5, linestyle="dotted")
    plt.savefig(synplot)
