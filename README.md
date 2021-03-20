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

The `.chrom` file format has 5 tab-delimited fields. Each line details the location of a protein on a scaffold. The fields are:
  - "protein_header scaffold_header strand start stop"

The requirements for each field are:
1. "protein_header" - the string here must match the header of a protein in the protein fasta.
2. "scaffold_header" - the string here must match the header of a sequence in the genome assembly fasta.
3. "strand" - must be "+" or "-".
4. "start" - the position, in basepair coordinates, where the CDS of the protein starts. Can often be found in a GFF3 or GTF.
5. "stop" - same as #4, but the stop position.

For example, the following text details four proteins that map to two scaffolds. Two of the proteins are on the negative strand. The first protein, `BFGX8636T1`, has its start codon from the first position of scaffold 1, and the last codon ends at base 1246.

```
BFGX8636T1      sca1    +       1       1246
BFGX0001T1      sca1    -       2059    2719
BFGX0002T1      sca2    +       6491    12359
BFGX0003T1      sca2    -       12899   18848
```

### How to genereate a `.chrom` file

A `.chrom` file can usually easily be generated from a genome annotation, such as a [`GFF3`](https://m.ensembl.org/info/website/upload/gff3.html) or [`GTF/GFF2`](https://uswest.ensembl.org/info/website/upload/gff.html) file.

