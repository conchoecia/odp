#!/usr/bin/env python

"""
This script tests the rbh_to_distance_gbgz script and how much ram it uses
"""

#import os
#import sys
# get the path of the current python script
#scriptpath = os.path.dirname(os.path.realpath(__file__))
#sys.path.insert(1, dependencies_path)
import PhyloTreeUMAP as PTU

rbh_file = "/lisc/scratch/molevo/dts/manifold/UMAP_snakemake_small/small_rbh/BCnSSimakov2022_Zalophuscalifornianus-9704-GCA009762305.2_xy_reciprocal_best_hits.plotted.rbh"
gbgz = "test.gb.gz"
ALGname = "BCnSSimakov2022"

PTU.rbh_to_distance_gbgz(rbh_file, gbgz, ALGname)
