#!/usr/bin/env python

"""This is for testing the plot of the results of annotating all the genomes."""

import os
import sys
# get the path of the current python script
scriptpath = os.path.dirname(os.path.realpath(__file__))
dependencies_path = os.path.join(scriptpath, "../dependencies/fasta-parser")
sys.path.insert(1, dependencies_path)
import AnnotateSampleDf as asd

df_filepath = "/lisc/scratch/molevo/dts/manifold/df_annotate/dfannotate/allsamples.neighbors_150.mind_0.9.missing_large.supplemented.df"
title       = "Chordates"
taxids_to_include = [7711]
taxids_to_exclude = []
pdfout            = "test_subclade_highlight.pdf"

asd.plot_UMAP_highlight_subclade(df_filepath,
                                 title,
                                 taxids_to_include,
                                 taxids_to_exclude,
                                 pdfout)
