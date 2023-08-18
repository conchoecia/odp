#!/usr/bin/env python3
"""
This program parses a NCBI GFF annotation and generates a .chrom file.
see https://github.com/conchoecia/odp for the specification.
"""
#NW_011887297.1	RefSeq	CDS	1566678	1566739	.	+	2	ID=cds-XP_004348322.1;Parent=rna-XM_004348272.1;Dbxref=GeneID:14898863,Genbank:XP_004348322.1;Name=XP_004348322.1;gbkey=CDS;locus_tag=CAOG_04494;product=hypothetical protein;protein_id=XP_004348322.1

import argparse
import datetime
import gzip
import os
import pandas as pd
import sys

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

# use the argparse library to get the gff filepath from the user
def parse_args():
    parser = argparse.ArgumentParser(description='Parse a NCBI GFF annotation and generate a .chrom file.')
    parser.add_argument('-g', '--gff',     type=str, required = True, help='Path to the GFF file.')
    # add a flag for the genome fasta
    parser.add_argument('-f', '--fasta',   type=str, required = True, help='Path to the genome fasta file.')
    # add a flag for the protein fasta
    parser.add_argument('-p', '--protein', type=str, required = True, help='Path to the protein fasta file.')
    # add a flag for the output prefix
    parser.add_argument('-o', '--outprefix',  type=str, required = True, help='Prefix for the output files.')
    args = parser.parse_args()

    # ensure that the files exists
    for thisfile in [args.gff, args.fasta, args.protein]:
        if not(os.path.isfile(thisfile)):
            parser.error("The file {} does not exist!".format(args.gff))

    return parser.parse_args()

def get_gff_handle(gff_path):
    """
    Reuse the code above to check if the file is gzipped or not
    """
    handle = None
    gzipped = False
    for thisend in [".gz", ".gzip", ".GZ", ".GZIP", ".gzipped", ".GZIPPED"]:
        if gff_path.endswith(thisend):
            gzipped = True

    if gzipped:
        handle = gzip.open(gff_path,'rt')
    else:
        handle = open(gff_path, "r")

    return handle

def get_fasta_headers_and_lengths(fastapath):
    """
    # Get the scaffold names from the fasta file and the scaffold length. 
    """
    scaf_to_len = {}
    # use SeqIO to get the headers and lengths
    # raise an error if the scaffold name appears more than once in the fasta file
    with open(fastapath, "r") as handle:
        for record in fasta.parse(handle):
            if record.id not in scaf_to_len:
                scaf_to_len[record.id] = len(record.seq)
            else:
                raise IOError("The sequence name {} appears more than once in the fasta file!".format(record.id))
    return scaf_to_len

def print_chrom_report(prots, scafnames_to_size, protnames_to_size, outprefix):
    """
    # Let's write the report first. We want to know how many proteins were in the gff file,
    #  how many proteins were in the headers passed in, and how many proteins were not in the gff file.
    #
    # In the report we will perform a similar analysis for the scaffolds. We need to know which scaffolds
    #  were not seen in the annotation.
    #
    # In the report we will also print the number of proteins that were on each scaffold.
    # We also will report the bias of forward or reverse strand.
    """
    report_path = "{}.report.txt".format(outprefix)

    # MAKE THE REPORT STRING FIRST
    # outstring is just called s for brevity
    s  = "#####  Gff3 to chrom conversion report.\n"
    s += "#####\n"
    s += "#####  Code that generated this script: https://github.com/conchoecia/odp/blob/main/scripts/NCBIgff2chrom.py\n"
    s += "#####  script author: dts (github @conchoecia)\n"
    s += "\n"
    s += "# Filename: {}\n".format(report_path)
    s += "# Date: {}\n".format(datetime.datetime.now())
    s += "#\n"
    protname_set = set(protnames_to_size.keys())
    prots_in_pep_but_not_gff = protname_set - set(prots.keys())
    s += "# PROTEINS:\n"
    s += "#   - There were {} proteins in the protein fasta file.\n".format(len(protnames_to_size))
    s += "#   - There were {} proteins in the gff file.\n".format(len(prots))
    s += "#   - Therefore, {} proteins were present in the protein fasta file but not in the gff file.\n".format(len(prots_in_pep_but_not_gff))
    s += "#   - Currently, this program does not run to completion if there are proteins in the gff file that are not in the protein fasta file.\n"
    s += "#\n"
    num_scafs_genome = len(scafnames_to_size)
    scafname_set = set(scafnames_to_size.keys())
    scafs_in_gff_file   = set([prots[x]["scaf"] for x in prots])
    scafs_in_genome_but_not_gff = scafname_set - scafs_in_gff_file
    s += "# SCAFFOLDS:\n"
    s += "#   - There were {} scaffolds in the genome fasta file.\n".format(len(scafname_set))
    s += "#   - There were {} scaffolds in the gff file.\n".format(len(scafs_in_gff_file))
    s += "#   - Therefore, {} scaffolds were present in the genome fasta file but not in the gff file.\n".format(len(scafs_in_genome_but_not_gff))
    s += "#   - Currently, this program does not run to completion if there are scaffolds in the gff file that are not in the genome fasta file.\n"
    s += "#\n"
    # now we need to calculate how many proteins were on each scaffold
    scaf2protcount = {}
    scaf2forwardcount = {}
    scaf2reversecount = {}
    for thisprot in prots:
        thisscaf = prots[thisprot]["scaf"]
        if thisscaf not in scaf2protcount:
            scaf2protcount[thisscaf] = 0
        scaf2protcount[thisscaf] += 1
        if prots[thisprot]["strand"] == "+":
            if thisscaf not in scaf2forwardcount:
                scaf2forwardcount[thisscaf] = 0
            scaf2forwardcount[thisscaf] += 1
        elif prots[thisprot]["strand"] == "-":
            if thisscaf not in scaf2reversecount:
                scaf2reversecount[thisscaf] = 0
            scaf2reversecount[thisscaf] += 1
    # turn scaf2protcount into a pandas df. The scaffold names will be a column, not an index
    df = pd.DataFrame(scaf2protcount.items(), columns=["scaf", "num_proteins"])
    df["scaflen"] = df["scaf"].map(scafnames_to_size)
    df["percent_of_proteins"] = 100 * (df["num_proteins"] / len(prots))
    df["forward_count"] = df["scaf"].map(scaf2forwardcount).fillna(0)
    df["reverse_count"] = df["scaf"].map(scaf2reversecount).fillna(0)
    df["FR_ratio"]      = df["forward_count"] / df["reverse_count"]
    # set types back to ints
    for thiscol in ["forward_count", "reverse_count"]:
        df[thiscol] = df[thiscol].astype(int)
    # sort the df by the number of proteins, descending 
    df = df.sort_values(by=["num_proteins"], ascending=False)
    s += "# PROTEINS PER SCAFFOLD:\n"
    s += df.to_string(index=False)
    s += "\n#\n"

    if os.path.exists(report_path):
        raise IOError("The report file {} already exists! This program will not overwrite existing files.".format(report_path))
    else:
        basedir = os.path.dirname(report_path)
        # safely make the directory for the report if it does not yet exist
        create_directories_recursive_notouch(basedir)
    
        outhandle = open(report_path, "w")
        print(s, file = outhandle)
        outhandle.close()
        # DONE WITH THE REPORT
        # give a safe return value
        return 0

def format_prots_to_df(prots):
    """
    This formats a dictionary we used earlier to store the protein information into a pandas dataframe for sorting and easy printing.

    The input format is like this:
        prots[pid] = {"scaf": scaf, "strand": strand,
                      "start": start, "stop": stop}

    We need to get it into a df with the same layout as a .chrom file
        protein ID, scaffold, strand, start, stop
    """
    entries = [{"pid": x, "scaf": prots[x]["scaf"], "strand": prots[x]["strand"], "start": prots[x]["start"], "stop": prots[x]["stop"]} for x in prots]
    df = pd.DataFrame(entries)
    del entries
    df = df.sort_values(by=["scaf", "start", "stop"])
    return df

def gff_to_chrom(gffhandle, genome_headers_to_size, protein_headers_to_size, outprefix):
    """
    A chrom file has the following format:
    protein_id	scaffold	strand	start	stop

    The protein_id must be present in the protein fasta file.
    The scaffold must be present in the genome fasta file.

    The input parameters are:
      - gffhandle: a handle to the gff file
      - genome_headers_to_size:  an iterable of the scaffold names in the genome fasta file
      - protein_headers_to_size: an iterable of the protein names in the protein fasta file
      - outprefix: the prefix for the output files. This script saves a .chrom and .pep file.
         N.B. - this script does not save over existing files. The program will throw an error
            if the output files already exist.
    
    The function saves the files:
      - .chrom file (safely saved to disk - no clobber)
      - .report.txt (safely saved to disk - no clobber)
    
    The function returns:
      - a pandas dataframe of the chrom file. This can be used to filter proteins later.
    """
    prots = {}

    # PARSE THE GFF FILE
    for line in gffhandle:
        line = line.strip()
        splitd=line.split("\t")
        if line and len(splitd) > 7 and splitd[2] == "CDS" and "protein_id=" in line:
            pid = [x for x in splitd[8].split(";") if x.startswith("protein_id=")][0].replace("protein_id=", "")
            scaf = splitd[0]
            # if this protein isn't in the protein fasta file, raise an error
            if pid not in protein_headers_to_size:
                raise IOError("The protein {} is not in the protein fasta file! Did you get the same annotation and protein file from NCBI?".format(pid))
            if scaf not in genome_headers_to_size:
                raise IOError("The scaffold {} is not in the genome fasta file! Did you get the same annotation and genome file from NCBI?".format(scaf))
            strand = splitd[6]
            start = int(splitd[3])
            stop = int(splitd[3])
            if pid not in prots:
                prots[pid] = {"scaf": scaf, "strand": strand,
                              "start": start, "stop": stop}
            else:
                if start < prots[pid]["start"]:
                    prots[pid]["start"] = start
                if stop > prots[pid]["stop"]:
                    prots[pid]["stop"] = stop
    # we don't need the gff file anymore
    gffhandle.close()

    # print the report
    print_chrom_report(prots, genome_headers_to_size, protein_headers_to_size, outprefix)

    df = format_prots_to_df(prots)

    # now print the chrom file
    chrom_path = "{}.chrom".format(outprefix) 
    if os.path.exists(chrom_path):
        raise IOError("The chrom file {} already exists! This program will not overwrite existing files.".format(report_path))
    else:
        basedir = os.path.dirname(chrom_path)
        # safely make the directory for the report if it does not yet exist
        create_directories_recursive_notouch(basedir)
    
        # save
        df = format_prots_to_df(prots)
        df.to_csv(chrom_path, sep="\t", index=False)
        # DONE WITH THE REPORT
        # give a safe return value
        return 0

def main():
    # first parse the args
    args = parse_args()

    # get the scaffold names from the fasta file 
    scafnames = get_fasta_headers_and_lengths(args.fasta)

    # get the protein names from the protein fasta file
    prots = get_fasta_headers_and_lengths(args.protein)

    # now get a handle for the gff file
    handle = get_gff_handle(args.gff)

    # make and print a report. The report saves to <outprefix>.report.txt
    # This also makes a protein file saved to <outprefix>.pep
    #  and it makes a chrom file saved to <outprefix>.chrom
    chromdf = gff_to_chrom(handle, scafnames, prots, args.outprefix)

if __name__ == '__main__':
    main()