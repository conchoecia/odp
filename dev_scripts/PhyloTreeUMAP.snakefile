"""
This is a snakefile that is used to parallelize constructing a chromosomal linkage UMAP.
 For this UMAP, each dot is a genome or an inferred pseudogenome at one node.

I am writing a snakefile to do this, because it takes too long in a for loop.
"""

configfile: "config.yaml"
# First we define all of the RBH files
# get the rbh files in the directory
rbh_files = list(sorted([os.path.join(args.directory, f)
             for f in os.listdir(args.directory)
             if f.endswith('.rbh')], reverse = True))
rbh_files = rbh_files[:100]


rule all:
