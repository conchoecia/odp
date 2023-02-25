#!/usr/bin/env python

import pandas as pd

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
    finaldf.to_csv(outfile, sep="\t", index = False, header = False)

reciprocal_best_permissive_blastp_or_diamond_blastp("COW_against_HCA.blastp", "HCA_against_COW.blastp", "HCA_against_COW.recipbest.blastp")
