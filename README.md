<img src="https://github.com/conchoecia/odp/blob/main/figs/dotplot_fig-01.jpg" width="250">

# odp - Oxford Dot Plots

This software creates dotplots of protein synteny in genome assemblies. This software was visually modelled on the dotplots found in [Simakov, Oleg, et al. "Deeply conserved synteny resolves early events in vertebrate evolution." Nature ecology & evolution 4.6 (2020): 820-830.](https://www.nature.com/articles/s41559-020-1156-z)

This software works by:
1. Finding reciprocal-best protein matches using diamond blastp.
2. Calculating [Da and Db](https://www.nature.com/articles/s41559-020-1156-z) for each genome.
3. Plotting the genome assembly, reciprocal best protein hits, and Da/Db using matplotlib

## Requirements

This software is probably best run in a unix environment. All development was on Ubuntu 18.04.4 LTS.

### python

Your active python environment must be python 3. This software is implemented in [`snakemake`](https://snakemake.readthedocs.io/en/stable/). Specific python packages within the pipeline that must be available in your path are:

```
snakemake
matplotlib
pandas
numpy
```
If you have `conda` I recommend `conda install snakemake matplotlib pandas numpy` if you are not sure if you have the required packages.

### other software

Direct calls to these programs must also be available in your environment. Future versions of `odp` may bundle these software packages directly to avoid these requirements.

- [diamond](https://github.com/bbuchfink/diamond)
- awk
- [pauvre](https://github.com/conchoecia/pauvre) - I may remove this requirement in future releases 
- [bioawk](https://github.com/lh3/bioawk) - I may remove this requirement in future releases


## Installation

```
git clone https://github.com/conchoecia/odp.git
```

## Input files

Three files are needed for each genome for which you would like to be plotted.
1. A genome assembly in [`.fasta` format](https://en.wikipedia.org/wiki/FASTA_format).
2. Protein sequences in [`.fasta` format](https://en.wikipedia.org/wiki/FASTA_format).
3. A `.chrom` file, which details where the proteins are located in the genome.

### `.chrom` specification

The `.chrom` file format has 5 fields. Each line details the location of a protein in the genome assem
