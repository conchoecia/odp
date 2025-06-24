#!/usr/bin/env python

"""This is for testing the plot of the results of annotating all the genomes."""

import os
import sys
# get the path of the current python script
scriptpath = os.path.dirname(os.path.realpath(__file__))
dependencies_path = os.path.join(scriptpath, "../dependencies/fasta-parser")
sys.path.insert(1, dependencies_path)
import AnnotateSampleDf as asd

alg_rbh = "/lisc/scratch/molevo/dts/odp/LG_db/BCnSSimakov2022/BCnSSimakov2022.rbh"
df      = "/lisc/scratch/molevo/dts/manifold/df_annotate/dfannotate/allsamples.neighbors_20.mind_0.9.missing_small.supplemented.df"
outpdf  = "annotatenew.pdf"
ALGname = "BCnSSimakov2022"

asd.bin_and_plot_decay(alg_rbh, df, outpdf, ALGname, 5)