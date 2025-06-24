#!/usr/bin/env python

"""
This is a test script for the function: umapdf_reimbedding_bokeh_plot_one_species
"""

from odol_annotate_blast import (umapdf_reimbedding_bokeh_plot_one_species,
                                 tsvgz_calcUMAP)

## Test 1, for deubgging and testing this function
## I think I am done with this testing
#filepath = "/lisc/scratch/molevo/dts/manifold/UMAP_blast/UMAP_blast_results/ALG_reimbedding/one_analysis_one_species_one_query/Coleoidea_6606_without_None.Octopusvulgaris-6645-GCA951406725.2.HOXA.B2.neighbors_15.mind_0.1.missing_large.subchrom.df"
#umapdf_reimbedding_bokeh_plot_one_species(filepath,
#                                          "Octopus vulgaris HOXA in ALG B2, nn=15, md=0.1",
#                                          "Octopus_vulgaris_HOXA_in_ALG_B2_nn15_md01.html")


# Test 2, testing the ALG_reimbedding function
# This takes about 18 seconds to run on its own.
# Octopodiformes!!
analysis   = "Octopodiformes_215451_without_None"
sampledf   = "/lisc/scratch/molevo/dts/manifold/UMAP_snakemake/GTUMAP/subchrom/phylogenetic/missing_large/Octopodiformes_215451_without_None.method_phylogenetic.missing_large.sampledf.tsv"
algrbhfile = "/lisc/scratch/molevo/dts/odp/LG_db/BCnSSimakov2022/BCnSSimakov2022.rbh"
tsvgz      = "/lisc/scratch/molevo/dts/manifold/UMAP_blast/UMAP_blast_results/filt_coo/Octopodiformes_215451_without_None_B2.matrix.filt.tsv.gz"
n          = 15
m          = 0.1
UMAPdf     = "Octopodiformes_215451_without_None.neighbors_15.mind_0.1.missing_large.subchrom.df"
UMAPbokeh  = "Octopodiformes_215451_without_None.neighbors_15.mind_0.1.missing_large.subchrom.bokeh.html"

tsvgz_calcUMAP(analysis, sampledf,
               algrbhfile, tsvgz,
               "large", n, m,
               UMAPdf)

# Test 2, testing the ALG_reimbedding function
# This takes about 18 seconds to run on its own.
# Octopodiformes!!
analysis   = "Octopodiformes_215451_without_None"
sampledf   = "/lisc/scratch/molevo/dts/manifold/UMAP_snakemake/GTUMAP/subchrom/phylogenetic/missing_large/Octopodiformes_215451_without_None.method_phylogenetic.missing_large.sampledf.tsv"
algrbhfile = "/lisc/scratch/molevo/dts/odp/LG_db/BCnSSimakov2022/BCnSSimakov2022.rbh"
tsvgz      = "/lisc/scratch/molevo/dts/manifold/UMAP_blast/UMAP_blast_results/filt_coo/Octopodiformes_215451_without_None_B2.matrix.filt.tsv.gz"
n          = 15
m          = 0.1
UMAPdf     = "Octopodiformes_215451_without_None.neighbors_15.mind_0.1.missing_large.subchrom.df"
UMAPbokeh  = "Octopodiformes_215451_without_None.neighbors_15.mind_0.1.missing_large.subchrom.bokeh.html"

tsvgz_calcUMAP(analysis, sampledf,
               algrbhfile, tsvgz,
               "large", n, m,
               UMAPdf)
