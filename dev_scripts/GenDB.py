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
from io import StringIO
from random import randint
import shutil
import subprocess
import sys
import time
import warnings

# use the fasta package included with this repo
snakefile_path = os.path.dirname(os.path.realpath(__file__))
dependencies_path = os.path.join(snakefile_path, "../dependencies/fasta-parser")
sys.path.insert(1, dependencies_path)
import fasta

def create_directories_recursive_notouch(path):
    """
    Unlike os.makedirs, this function will not touch a directory if it already exists.
    This is useful for snakemake because it will not re-run a rule if the output already exists.
    """
    parts = os.path.normpath(path).split(os.path.sep)
    current_path = ""

    for part in parts:
        current_path = os.path.join(current_path, part)
        if not os.path.exists(current_path):
            os.mkdir(current_path)

def gzip_get_time(basepairs) -> int:
    """
    This function returns how many minutes of compute time a gzip job on a fasta file will take.
    This was inferred from testing around 2400 genomes of lengths between 100 Mbp to 40 Gbp.
    For every Gbp of bases, it takes about 7 minutes of gzip time. We add an additional 2 minutes of overhead.

    The equation in the form y=mx+b is:
       time = 7.0 * (basepairs/1e9) + 2.0

    Returns an int, because SLURM only takes ints for runtime.
    """
    return int(7.0 * (basepairs/1e9) + 2.0)

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
        #if chosen_file is None:
        #    raise IOError("You provided a directory for us to find the lists of genomes {}, but")
        # Use the chosen file, and add the entries to the config file so we can download them or not
        df = pd.read_csv(chosen_file, sep="\t", index_col=False)
        print(df)
        # strip leading and trailing whitespace from the column names because pandas can screw up sometimes
        df.columns = df.columns.str.strip()

        config["assemAnn"]            = df["Assembly Accession"].tolist()
        # make a dict where assemAnn is the key and the value is the "Assembly Stats Total Ungapped Length" column
        config["assemAnn_to_scaflen"] = dict(zip(df["Assembly Accession"].tolist(), df["Assembly Stats Total Ungapped Length"].tolist()))

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

def download_assembly_scaffold_df(assembly_accession, datasetsEx, chrscale = True) -> pd.DataFrame:
    """
    Handles downloading the assembly scaffold df from NCBI.
    Uses the NCBI datasets executable tool to do so.

    Input:
        - assembly_accession: the assembly accession number of the genome to download
        - datasetsEx: path to the NCBI datasets executable
    Output:
        - df: a pandas dataframe containing the assembly scaffold df

    df output example:
                assembly_accession     assembly_unit assigned_molecule_location_type chr_name    gc_count  gc_percent  genbank_accession    length                role
        0      GCA_026151205.1  Primary Assembly                      Chromosome        1  10152078.0        38.5         CM047403.1  26206493  assembled-molecule
        1      GCA_026151205.1  Primary Assembly                      Chromosome        2   8914563.0        38.5         CM047404.1  23028209  assembled-molecule
        2      GCA_026151205.1  Primary Assembly                      Chromosome        3   8250754.0        39.0         CM047405.1  21284423  assembled-molecule
        3      GCA_026151205.1  Primary Assembly                      Chromosome        4   6284100.0        38.5         CM047406.1  16326010  assembled-molecule
        4      GCA_026151205.1  Primary Assembly                      Chromosome        5   9163041.0        39.0         CM047407.1  23637599  assembled-molecule
        ..                 ...               ...                             ...      ...         ...         ...                ...       ...                 ...
        342    GCA_026151205.1  Primary Assembly                      Chromosome       Un         NaN         NaN  JAIOUN010000343.1      3450   unplaced-scaffold
        343    GCA_026151205.1  Primary Assembly                      Chromosome       Un         NaN         NaN  JAIOUN010000344.1      3041   unplaced-scaffold
        344    GCA_026151205.1  Primary Assembly                      Chromosome       Un         NaN         NaN  JAIOUN010000345.1      2526   unplaced-scaffold
        345    GCA_026151205.1  Primary Assembly                      Chromosome       Un         NaN         NaN  JAIOUN010000346.1      1772   unplaced-scaffold
        346    GCA_026151205.1       non-nuclear                   Mitochondrion       MT      2109.0        20.0         CM047416.1     10476  assembled-molecule
    """
    # wait a random amount of time up to 120 seconds to space out requests to the NCBI servers.
    sleeptime = randint(1,120)
    print("About to download the assembly df for {}.\nSleeping for {} seconds to avoid overloading the NCBI servers.".format(assembly_accession,sleeptime), file = sys.stderr)
    time.sleep(sleeptime)

    # first use the datasets executable to get a list of scaffolds to download
    # launch it as a process from bash and collect the output to sys.stdout
    # ./datasets summary genome accession GCF_000001405.40 --report sequence --as-json-lines
    cmd = [datasetsEx, "summary", "genome", "accession", assembly_accession, "--report", "sequence", "--as-json-lines"]
    p = ""
    if chrscale:
        # if we are downloading the chromosome-scale scaffolds, then we can pre-filter the non-chromosome-scale scaffolds
        # add this to the command: | grep -v 'unplaced-scaffold'
        cmd += ["|", "grep", "-v", "'unplaced-scaffold'", "|", "grep", "-v", "'chr_name\":\"Un'"]
        cmd = " ".join(cmd)
        print("Trying the following command: {}".format(cmd), file = sys.stderr)
        p = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    else:
        print("Trying the following command: {}".format(" ".join(cmd)), file = sys.stderr)
        p = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    out, err = p.communicate()
    sys.stderr.write(err)
    if p.returncode != 0:
        raise ValueError("Unzip for {} returned a non-zero exit code. Something went wrong. {}".format(assembly_accession, err))
    # the output is json lines. Convert it to a pandas dataframe
    json_buffer = StringIO(out)
    df = pd.read_json(json_buffer, lines=True)
    # return an error if the dataframe is empty
    if df.empty:
        raise ValueError("The dataframe is empty. This is probably because the assembly accession number you provided is not valid. {}".format(assembly_accession))
    # return a dataframe if there are any values in the "assembly_accession" column that are not the same as the one we provided
    if not df["assembly_accession"].unique()[0] == assembly_accession:
        raise ValueError("The assembly accession number you provided is not valid. {}".format(assembly_accession))
    # return the dataframe
    return df

def filter_scaffold_df_keep_chrs(df) -> pd.DataFrame:
    """
    Takes in an NCBI dataset of the scaffold entries and filters it down to a df that only contains the chromosome-scale scaffolds.
    The input dataframe looks like the one detailed in the docstring of download_assembly_scaffold_df(assembly_accession, datasetsEx)
    The strategy to get the chromosome-scale scaffolds is mostly to remove the things that we know aren't chromosome-scale.

    Input:
        - scaffold_df: a pandas dataframe containing the assembly scaffold df
    Output:
        - df: a pandas dataframe containing only the chromosome-scale scaffolds

    Notes:
        - DO NOT filter on "Chromosome" in the "assigned_molecule_location_type" column. This will remove linkage groups, which are also chromosome-scale scaffolds.
          Instead remove "Mitochondrion" entries from the "assigned_molecule_location_type" column.
    """
    startdf = df.copy()
    # get the most common assembly accession number in the df
    assembly_accession = df["assembly_accession"].mode()[0]
    start_length = len(df)

    # Remove mitochondria
    if "assigned_molecule_location_type" in df.columns:
        df = df[df["assigned_molecule_location_type"] != "Mitochondrion"]
    if "role" in df.columns:
        # remove things called "unplaced-scaffold" - these are non-chromosome-scale scaffolds for which the chromosome is uncertain
        df = df[df["role"] != "unplaced-scaffold"]
        # remove things called "alt-scaffold" - these are not the main scaffolds for the annotation
        df = df[df["role"] != "alt-scaffold"]
        # remove things that are called "unlocalized-scaffold", the chromosome to which these belond is known, but the location on the chromosome is not known
        df = df[df["role"] != "unlocalized-scaffold"]
    if "chr_name" in df.columns:
        # remove chr_name Un
        df = df[df["chr_name"] != "Un"]
    # raise an error if the dataframe is empty - we expect there to be at least one chromosome-scale scaffold
    if df.empty:
        print("This was the startdf: \n{}".format(startdf), file = sys.stderr)
        raise ValueError("The dataframe is empty after filtering for chromosome-scale scaffolds. The start length of the df was {}. The assembly_accession is {}".format(start_length, assembly_accession))
    return df

def download_chr_scale_genome_from_df_WORKING(chr_df, datasetsEx, output_dir):
    """
    NOTE: THIS FUNCTION IS ONLY TO BE USED UNTIL https://github.com/ncbi/datasets/issues/298 IS FIXED
    This function takes in a dataframe of known chromosome-scale scaffolds for a specific assembly accession number,
      and downloads the genome from NCBI.

    Special checks that are made along the way to make sure we get exactly what is needed:
      - We check that we only add each needed scaffold exactly once. This avoids a problem of using a bash-based solution
        where we cannot finely control which fasta files are concatenated together.
      - We check that each fasta entry is exactly the length that is expected. This will help find cases where the
        download or decompression did not work for some reason.

    Other features:
      - This function will delete all of the unneeded intermediate files.

    Inputs:
        - chr_df:      a pandas dataframe containing only the chromosome-scale scaffolds
        - datasetsEx:  path to the NCBI datasets executable
        - output_dir:  The directory to which additional files will be saved.
                         Everything will be prefixed with the assembly_accession value.
    Outputs:
        - Outputs 0 if the function completes successfully.
        - Saves a fasta file to the output_dir that contains all of the chromosome-scale scaffolds.
          The file will be saved at {output_dir}/{assembly_accession}.chr.fasta
    """
    # get the most common value of the assembly_accession column
    assembly_accession = chr_df["assembly_accession"].mode()[0]

    # wait a random amount of time up to 120 seconds to space out requests to the NCBI servers.
    sleeptime = randint(1,120)
    print("About to download the genome assembly for {}.\nSleeping for {} seconds to avoid overloading the NCBI servers.".format(assembly_accession, sleeptime), file = sys.stderr)
    time.sleep(sleeptime)

    # First, we need to check that the output directory exists. If it doesn't, then we need to create it.
    create_directories_recursive_notouch(output_dir)

    # DOWNLOAD THE SCAFFOLDS
    # We download the entire genome, and then filter it down to just the chromosome-scale scaffolds.
    # use the assembly accession number to name the output zip file
    outzip = os.path.join(output_dir, "{}.zip".format(assembly_accession))
    # set up the download process using the datasets executable
    cmd = [datasetsEx, "download", "genome", "accession", assembly_accession, "--filename", outzip]
    print("the command is {}".format(" ".join(cmd)), file = sys.stderr)
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text = True)
    # wait for the process to finish
    out, err = p.communicate()
    sys.stdout.write(out)
    sys.stderr.write(err)
    if p.returncode != 0:
        raise ValueError("The datasets executable for {} returned a non-zero exit code. Something went wrong. {}".format(assembly_accession, err))
    # check that the output zip file exists
    if not os.path.exists(outzip):
        raise ValueError("The output zip file does not exist. Something went wrong with the download. {}".format(outzip))

    # UNZIP THE SCAFFOLDS
    # unzip the file, forcing it to overwrite any existing files, excluding the file called README.md
    unzip_dir = os.path.join(output_dir, assembly_accession)
    cmd = ["unzip", "-o", outzip, "-d", unzip_dir]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text = True)
    # wait for the process to finish
    out, err = p.communicate()
    sys.stdout.write(out)
    sys.stderr.write(err)
    if p.returncode != 0:
        raise ValueError("Unzip for {} returned a non-zero exit code. Something went wrong. {}".format(assembly_accession, err))
    # check that the output directory exists
    if not os.path.isdir(unzip_dir):
        raise ValueError("The output directory does not exist. Something went wrong with the unzip. {}".format(output_dir))

    # CHECK THAT THE MAIN ASSEMBLY FILE IS THERE. It will be called GCF_905220415.1_{something}_genomic.fna
    found_file = False
    target_dir = os.path.join(unzip_dir, "ncbi_dataset/data/{}".format(assembly_accession))
    target_file = ""
    for thisfile in os.listdir(target_dir):
        if thisfile.endswith("_genomic.fna"):
            found_file = True
            target_file = os.path.join(target_dir, thisfile)
            break
    if not found_file:
        # The file does not exist. Something went wrong with the download of one of the chromosomes
        raise IOError("For Assembly Accession {}, the file {}/{}_something_genomic.fna does not exist. Something went wrong with the download of this chromosome.".format(assembly_accession, target_dir, thisfile))
    print(target_file)

    # CHECK THAT THE SCAFFOLDS ARE THE CORRECT LENGTH AND SAVE TO NEW FILE
    # To keep track of which scaffolds we need to include in the output file, use a copy of the chr_df and delete rows as we go.
    outfile = os.path.join(output_dir, "{}.chr.fasta".format(assembly_accession))
    outhandle = open(outfile, "w")
    scafs_to_parse_df = chr_df.copy()

    col_to_check = ""
    # go through the fasta file of the assembly, and look for whether the earliest scaffolds (chrs) are genbank_accession or refseq_accession
    # Based on this result, we will know which column to check
    scaf_to_col = {}
    for thiscol in ["genbank_accession", "refseq_accession"]:
        if thiscol in scafs_to_parse_df.columns:
            for thisscaf in scafs_to_parse_df[thiscol].tolist():
                scaf_to_col[thisscaf] = thiscol
    for record in fasta.parse(target_file):
        if record.id in scaf_to_col:
            col_to_check = scaf_to_col[record.id]
            break
        else:
            # This scaffold is not in the dataframe. We shouldn't be here, but maybe one unrequested scaffold was downloaded by the NCBI datasets tool.
            # One instance in which this happens often is when we tell datasets to download specific chromosomes, but it also downloads separate fasta files for the unplaced scaffolds.
            pass
    if col_to_check == "":
        raise ValueError("We didn't find any scaffolds in the fasta file that were in the dataframe. Something went wrong. Assembly accession: {}".format(assembly_accession))

    # Iterate through the files in the directory and:
    #  - go through each entry in the fasta file
    #    - check if the scaffold is in scafs_to_parse_df
    #    - if it is, then check that the length is correct
    #    - print it out to the outhandle
    #    - remove the corresponding row from scafs_to_parse_df
    for record in fasta.parse(target_file):
        # get the corresponding row from scafs_to_parse_df. could be either "genbank_accession" or "refseq_accession"
        if record.id in scafs_to_parse_df[col_to_check].tolist():
            # there is a row in the dataframe that corresponds to this scaffold
            # Get all the indices that correspond to this scaffold, using either the "genbank_accession" or "refseq_accession" column
            matching_indices = scafs_to_parse_df.index[scafs_to_parse_df[col_to_check] == record.id].tolist()
            if len(matching_indices) > 1:
                # There should only be one row that matches this scaffold. If there are more, then something went wrong.
                raise ValueError("There are multiple rows in the dataframe that match this scaffold. Something went wrong (Presumably with NCBI datasets, as this dataframe was generated with 'datasets summary genome accession ACC --report sequence --as-json-lines'). {}".format(record.id))
            # If we get here, then there is only one row that matches this scaffold. Get the index of that row.
            rowix = matching_indices[0]
            # Get the expected length of this scaffold, based on this row's "length" column
            expected_length = scafs_to_parse_df.loc[rowix, "length"]
            # After passing these checks, safe to print to the output file
            print(record.format(wrap=80), file=outhandle, end="")
            # remove the row from scafs_to_parse_df, reset the index
            scafs_to_parse_df = scafs_to_parse_df.drop(rowix).reset_index(drop=True)
        else:
            # This scaffold is not in the dataframe. We shouldn't be here, but maybe one unrequested scaffold was downloaded by the NCBI datasets tool.
            # One instance in which this happens often is when we tell datasets to download specific chromosomes, but it also downloads separate fasta files for the unplaced scaffolds.
            pass
    # close the output file
    outhandle.close()
    # check if there are any rows left in scafs_to_parse_df. If there are, then we didn't see one chromosome's sequence for some reason.
    if not scafs_to_parse_df.empty:
        # remove the output fasta file, since it is broken
        os.remove(outfile)
        # remove the unzipped directory and the zip file
        os.remove(outzip)
        shutil.rmtree(unzip_dir)
        raise ValueError("There are still rows in scafs_to_parse_df. This means that we didn't see one of the chromosomes. Something went wrong. {}".format(scafs_to_parse_df))
    # Just to be completely sure that we don't have any doubled scaffolds appearing in the output fasta file, read through it again with the fasta package
    seen_once = []
    for record in fasta.parse(outfile):
        if record.id in seen_once:
            raise ValueError("The scaffold {} appears twice in the output fasta file. Something went wrong.".format(record.id))
        else:
            seen_once.append(record.id)
    # Another excessive check, just make sure that the length of seen_once is the same length as the number of rows in the chr_df
    if not len(seen_once) == len(chr_df):
        raise ValueError("The number of scaffolds in the output fasta file is not the same as the number of scaffolds in the input dataframe. Something went wrong. Assembly acession: {}".format(assembly_accession))
    # If we get here, then everything probably worked as we intended. Remove the unneeded files.
    os.remove(outzip)
    shutil.rmtree(unzip_dir)
    return 0


def download_chr_scale_genome_from_df(chr_df, datasetsEx, output_dir):
    """
    NOTE: THIS WILL ONLY WORK WHEN https://github.com/ncbi/datasets/issues/298 IS FIXED
    This function takes in a dataframe of known chromosome-scale scaffolds for a specific assembly accession number,
      and downloads the genome from NCBI.

    Special checks that are made along the way to make sure we get exactly what is needed:
      - We check that we only add each needed scaffold exactly once. This avoids a problem of using a bash-based solution
        where we cannot finely control which fasta files are concatenated together.
      - We check that each fasta entry is exactly the length that is expected. This will help find cases where the
        download or decompression did not work for some reason.

    Other features:
      - This function will delete all of the unneeded intermediate files.

    Inputs:
        - chr_df:      a pandas dataframe containing only the chromosome-scale scaffolds
        - datasetsEx:  path to the NCBI datasets executable
        - output_dir:  The directory to which additional files will be saved.
                         Everything will be prefixed with the assembly_accession value.
    Outputs:
        - Outputs 0 if the function completes successfully.
        - Saves a fasta file to the output_dir that contains all of the chromosome-scale scaffolds.
          The file will be saved at {output_dir}/{assembly_accession}.chr.fasta
    """
    # get the most common value of the assembly_accession column
    assembly_accession = chr_df["assembly_accession"].mode()[0]
    # First, we need to check that the output directory exists. If it doesn't, then we need to create it.
    create_directories_recursive_notouch(output_dir)
    # cast the chr_name column to strings so we can work with them as strings later
    chr_df["chr_name"] = chr_df["chr_name"].astype(str)

    # DOWNLOAD THE SCAFFOLDS
    # Now we need to figure out which scaffolds we need to download. To do this, we can use the chromosome names and prompt the datasets tool to download them
    chrom_download_string = ",".join(chr_df["chr_name"].tolist())
    print("the chrom download string is {}".format(chrom_download_string), file = sys.stderr)
    # use the assembly accession number to name the output zip file
    outzip = os.path.join(output_dir, "{}.zip".format(assembly_accession))
    # set up the download process using the datasets executable
    cmd = [datasetsEx, "download", "genome", "accession", assembly_accession, "--chromosomes", chrom_download_string, "--filename", outzip]
    print("the command is {}".format(" ".join(cmd)), file = sys.stderr)
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text = True)
    # wait for the process to finish
    out, err = p.communicate()
    sys.stdout.write(out)
    sys.stderr.write(err)
    if p.returncode != 0:
        raise ValueError("The datasets executable for {} returned a non-zero exit code. Something went wrong. {}".format(assembly_accession, err))
    # check that the output zip file exists
    if not os.path.exists(outzip):
        raise ValueError("The output zip file does not exist. Something went wrong with the download. {}".format(outzip))

    # UNZIP THE SCAFFOLDS
    # unzip the file, forcing it to overwrite any existing files, excluding the file called README.md
    unzip_dir = os.path.join(output_dir, assembly_accession)
    cmd = ["unzip", "-o", outzip, "-d", unzip_dir]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text = True)
    # wait for the process to finish
    out, err = p.communicate()
    sys.stdout.write(out)
    sys.stderr.write(err)
    if p.returncode != 0:
        raise ValueError("Unzip for {} returned a non-zero exit code. Something went wrong. {}".format(assembly_accession, err))
    # check that the output directory exists
    if not os.path.isdir(unzip_dir):
        raise ValueError("The output directory does not exist. Something went wrong with the unzip. {}".format(output_dir))

    # CHECK THAT ALL THE SCAFFOLD FILES ARE THERE
    for thischrom in chr_df["chr_name"].tolist():
        thisfile = os.path.join(unzip_dir, "ncbi_dataset/data/{}/chr{}.fna".format(assembly_accession,thischrom))
        if not os.path.exists(thisfile):
            # The file does not exist. Something went wrong with the download of one of the chromosomes
            raise IOError("For Assembly Accession {}, the file {} does not exist. Something went wrong with the download of this chromosome.".format(assembly_accession, thisfile))

    # CHECK THAT THE SCAFFOLDS ARE THE CORRECT LENGTH AND SAVE TO NEW FILE
    # To keep track of which scaffolds we need to include in the output file, use a copy of the chr_df and delete rows as we go.
    outfile = os.path.join(output_dir, "{}.chr.fasta".format(assembly_accession))
    outhandle = open(outfile, "w")
    scafs_to_parse_df = chr_df.copy()
    # Iterate through the files in the directory and:
    #  - go through each entry in the fasta file
    #    - check if the scaffold is in scafs_to_parse_df
    #    - if it is, then check that the length is correct
    #    - print it out to the outhandle
    #    - remove the corresponding row from scafs_to_parse_df
    for thischrom in chr_df["chr_name"].tolist():
        thisfile = os.path.join(unzip_dir, "ncbi_dataset/data/{}/chr{}.fna".format(assembly_accession,thischrom))
        for record in fasta.parse(thisfile):
            # get the corresponding row from scafs_to_parse_df. could be either "genbank_accession" or "refseq_accession"
            if record.id in scafs_to_parse_df["genbank_accession"].tolist() + scafs_to_parse_df["refseq_accession"].tolist():
                # there is a row in the dataframe that corresponds to this scaffold
                # Get all the indices that correspond to this scaffold, using either the "genbank_accession" or "refseq_accession" column
                matching_indices = scafs_to_parse_df.index[scafs_to_parse_df["genbank_accession"] == record.id].tolist() + scafs_to_parse_df.index[scafs_to_parse_df["refseq_accession"] == record.id].tolist()
                if len(matching_indices) > 1:
                    # There should only be one row that matches this scaffold. If there are more, then something went wrong.
                    raise ValueError("There are multiple rows in the dataframe that match this scaffold. Something went wrong (Presumably with NCBI datasets, as this dataframe was generated with 'datasets summary genome accession ACC --report sequence --as-json-lines'). {}".format(record.id))
                # If we get here, then there is only one row that matches this scaffold. Get the index of that row.
                rowix = matching_indices[0]
                # Get the expected length of this scaffold, based on this row's "length" column
                expected_length = scafs_to_parse_df.loc[rowix, "length"]
                # if the expected length doesn't match the actual length, then something went wrong.
                # THIS IS COMMENTED OUT UNTIL https://github.com/ncbi/datasets/issues/301 IS FIXED
                #if not expected_length == len(record.seq):
                #    raise ValueError("The expected length of scaffold {} from assembly accession {} is {}, but the actual length is {}. Something went wrong.".format(record.id, assembly_accession, expected_length, len(record.seq)))
                # After passing these checks, safe to print to the output file
                print(record.format(wrap=80), file=outhandle, end="")
                # remove the row from scafs_to_parse_df, reset the index
                scafs_to_parse_df = scafs_to_parse_df.drop(rowix).reset_index(drop=True)
            else:
                # This scaffold is not in the dataframe. We shouldn't be here, but maybe one unrequested scaffold was downloaded by the NCBI datasets tool.
                # One instance in which this happens often is when we tell datasets to download specific chromosomes, but it also downloads separate fasta files for the unplaced scaffolds.
                pass
    # close the output file
    outhandle.close()
    # check if there are any rows left in scafs_to_parse_df. If there are, then we didn't see one chromosome's sequence for some reason.
    if not scafs_to_parse_df.empty:
        # remove the output fasta file, since it is broken
        os.remove(outfile)
        # remove the unzipped directory and the zip file
        os.remove(outzip)
        shutil.rmtree(unzip_dir)
        raise ValueError("There are still rows in scafs_to_parse_df. This means that we didn't see one of the chromosomes. Something went wrong. {}".format(scafs_to_parse_df))
    # Just to be completely sure that we don't have any doubled scaffolds appearing in the output fasta file, read through it again with the fasta package
    seen_once = []
    for record in fasta.parse(outfile):
        if record.id in seen_once:
            raise ValueError("The scaffold {} appears twice in the output fasta file. Something went wrong.".format(record.id))
        else:
            seen_once.append(record.id)
    # Another excessive check, just make sure that the length of seen_once is the same length as the number of rows in the chr_df
    if not len(seen_once) == len(chr_df):
        raise ValueError("The number of scaffolds in the output fasta file is not the same as the number of scaffolds in the input dataframe. Something went wrong. Assembly acession: {}".format(assembly_accession))
    # If we get here, then everything probably worked as we intended. Remove the unneeded files.
    os.remove(outzip)
    shutil.rmtree(unzip_dir)
    return 0

def download_unzip_genome(assembly_accession, output_dir, datasetsEx,
                          chrscale = True):
    """
    Given an assembly accession number, download the genome from NCBI and unzip it.
    If the assembly is chromosome-scale, then we:
      - do extra processing to only download the chromsome-scale scaffolds
      - unzip the genome
      - produce a single chromosome-scale genome file
      - check that the fasta file contains all of the chromsoome-scale scaffolds
      - check that all of the chromosome-scale scaffolds are the correct length

    Inputs:
      - assembly_accession:   The assembly accession number of the genome to download
      - output_dir:           The directory to which additional files will be saved.
                              Everything will be prefixed with the assembly_accession value.
      - final_fasta_filepath: The path to the final fasta file that we will produce
      - datasetsEx:           path to the NCBI datasets executable
      - dataformatEx:         path to the NCBI dataformat executable
      - chrscale:             boolean, whether the genome is chromosome-scale or not. Invokes different processing.

    Returns:
      - If a chromosome-scale genome, there will be a file called {assembly_accession}.chr.fasta in the output_dir.
      - If not a chromosome-scale genome, there will be a file called {assembly_accession}.all.fasta in the output_dir.
      - Returns 0 if the function completes successfully.
    """
    # First, strip the white space surrounding the assembly accession number
    assembly_accession = assembly_accession.strip()

    # save the current directory in case we need to go back to it
    current_dir = os.getcwd()
    if chrscale:
        scaffold_df = download_assembly_scaffold_df(assembly_accession, datasetsEx, chrscale = True)
        # we're going to save a file to a directory. First, make sure that the directory exists.
        create_directories_recursive_notouch(output_dir)
        # Originally, I was going to first check if the scaffold_df was already saved to output_dir.
        #  Then, I thought this would be a bad idea because the entry could change on NCBI's end, in which case it is best to download the scaffold info every time.
        scaffold_df_filepath = os.path.join(output_dir, assembly_accession + ".scaffold_df.all.tsv")
        # save the scaffold_df to a file
        scaffold_df.to_csv(scaffold_df_filepath, sep="\t", index=False)
        # Now we just get the chromosome-scale scaffolds
        chr_df = filter_scaffold_df_keep_chrs(scaffold_df)
        chr_df_filepath = os.path.join(output_dir, assembly_accession + ".scaffold_df.chr.tsv")
        chr_df.to_csv(chr_df_filepath, sep="\t", index=False)
        # Figure out which scaffolds we need to download. To do this, we can use the chromosome names and prompt the datasets tool to download them
        download_chr_scale_genome_from_df_WORKING(chr_df, datasetsEx, output_dir)
        # check that there is a file in the output
        if not os.path.exists(os.path.join(output_dir, assembly_accession + ".chr.fasta")):
            raise IOError("The output fasta file does not exist. Something went wrong with the download. {}".format(os.path.join(output_dir, assembly_accession + ".chr.fasta")))
    else:
        # The genome is not chromsome-scale, or we want to download every scaffold.
        pass
    return 0

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

    # It is fine if the files are empty, but there needs to at least be a file to look at
    cycledict = {"annotated_chr_files":      annotated_chr_files,
                    "annotated_nonchr_files":   annotated_nonchr_files,
                    "unannotated_chr_files":    unannotated_chr_files,
                    "unannotated_nonchr_files": unannotated_nonchr_files}
    for thiskey in cycledict:
        if len(cycledict[thiskey]) == 0:
            outmsg = "There are no accession TSV files for this filetype: {}".format(thiskey)
            outmsg = outmsg + "\nThe directory in which I looked was: {}".format(directory_path)
            outmsg = outmsg + "\n  - The files that we found in the directory are:"
            for thisfile in os.listdir(directory_path):
                outmsg = outmsg + "\n    - {}".format(thisfile)
            outmsg = outmsg + "\nThe names of the files that we need to find are found like so:"
            outmsg = outmsg + "\n    [f for f in files if f.startswith(\"annotated_genomes_chr\")      and f.endswith(\".tsv\") and contains_date(f)]"
            outmsg = outmsg + "\n    [f for f in files if f.startswith(\"annotated_genomes_nonchr\")   and f.endswith(\".tsv\") and contains_date(f)]"
            outmsg = outmsg + "\n    [f for f in files if f.startswith(\"unannotated_genomes_chr\")    and f.endswith(\".tsv\") and contains_date(f)]"
            outmsg = outmsg + "\n    [f for f in files if f.startswith(\"unannotated_genomes_nonchr\") and f.endswith(\".tsv\") and contains_date(f)]"
            outmsg = outmsg + "\nIt is fine if the files are empty, but there needs to at least be a file to look at."
            warnings.warn(outmsg) # using a warning for now, as it is not a fatal error. In case something is missing, this will just not run.

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