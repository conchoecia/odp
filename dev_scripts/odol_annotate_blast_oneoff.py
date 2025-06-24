#!/usr/bin/env python

# We need a few plots, just run them individually.

from odol_annotate_blast import umapdf_reimbedding_bokeh_plot_one_species


blastdf = "/lisc/scratch/molevo/dts/manifold/UMAP_blast/UMAP_blast_results/ALG_reimbedding/Vertebrata_7742_without_None/B2/Vertebrata_7742_without_None.neighbors_15.mind_0.1.missing_large.subchrom.df"
plot_title  = "vertebrate HoxA mind = 0.1, nneighbors = 15"
output_html = "Vertebrata_7742_without_None.neighbors_15.mind_0.1.missing_large.subchrom.df.bokeh.html"
umapdf_reimbedding_bokeh_plot_one_species(blastdf, plot_title, output_html, scalar = 1.5)

blastdf = "/lisc/scratch/molevo/dts/manifold/UMAP_blast/UMAP_blast_results/ALG_reimbedding/Echinodermata_7586_without_None/B2/Echinodermata_7586_without_None.neighbors_15.mind_0.1.missing_large.subchrom.df"
plot_title  = "Echinoderm HoxA mind = 0.1, nneighbors = 15"
output_html = "Echinodermata_7586_without_None.neighbors_15.mind_0.1.missing_large.subchrom.df.bokeh.html"
umapdf_reimbedding_bokeh_plot_one_species(blastdf, plot_title, output_html, scalar = 1.5)
