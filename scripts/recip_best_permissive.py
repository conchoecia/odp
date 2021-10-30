#!/usr/bin/env python

import pandas as pd

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
