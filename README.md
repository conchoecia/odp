<img src="https://github.com/conchoecia/odp/blob/main/figs/dotplot_fig-01.jpg" width="250">

# odp - Oxford Dot Plots

This software creates dotplots of protein synteny in genome assemblies. This software was visually modelled on the dotplots found in [Simakov, Oleg, et al. "Deeply conserved synteny resolves early events in vertebrate evolution." Nature ecology & evolution 4.6 (2020): 820-830.](https://www.nature.com/articles/s41559-020-1156-z)

This software works by:
1. Finding reciprocal-best protein matches using diamond blastp.
2. Calculating [Da and Db](https://www.nature.com/articles/s41559-020-1156-z) for each genome.
3. Plotting the genome assembly, reciprocal best protein hits, and Da/Db using matplotlib

## Example


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
4. "start" - the numerically least position, in basepair coordinates, of the CDS of the protein. Like the start coords in a `bed` or `GFF` file, not necessarily the position of the start codon. Can often be found in a GFF3 or GTF.
5. "stop" - same as #4, but the stop position.

For example, the following text details four proteins that map to two scaffolds. Two of the proteins are on the negative strand. The first protein, `BFGX8636T1`, has its start codon from the first position of scaffold 1, and the last codon ends at base 1246.

```
BFGX8636T1      sca1    +       1       1246
BFGX0001T1      sca1    -       2059    2719
BFGX0002T1      sca2    +       6491    12359
BFGX0003T1      sca2    -       12899   18848
```

### How to genereate a `.chrom` file

A `.chrom` file can usually easily be generated from a genome annotation, such as a [`GFF3`](https://m.ensembl.org/info/website/upload/gff3.html) or [`GTF/GFF2`](https://uswest.ensembl.org/info/website/upload/gff.html) file. If you are working with NCBI GFFs, CDS entries have a predictable format that enables us to compile all of the information required for a chrom file: `NC_000001.11	BestRefSeq	CDS	65565	65573	.	+	0	ID=cds-NP_001005484.2;Parent=rna-NM_001005484.2;Dbxref=CCDS:CCDS30547.1,Ensembl:ENSP00000493376.2,GeneID:79501,Genbank:NP_001005484.2,HGNC:HGNC:14825;Name=NP_001005484.2;gbkey=CDS;gene=OR4F5;product=olfactory receptor 4F5;protein_id=NP_001005484.2;tag=MANE Select`. There are many CDS lines per gene, so a special parsing program is required.

The program bundled with `odp`, [scripts/NCBIgff2chrom.py](scripts/NCBIgff2chrom.py), parses gzipped/uncompressed NCBI GFFs and gets the full protein span in the genome. Running [scripts/NCBIgff2chrom.py](scripts/NCBIgff2chrom.py) on the [human GFF from NCBI](https://www.ncbi.nlm.nih.gov/genome/?term=human) with the command `python NCBIgff2chrom.py GCF_000001405.39_GRCh38.p13_genomic.gff.gz` results in a legal `.chrom` file with all of the proteins from the annotation. This file can be easily filtered later on.

```
NP_001005484.2  NC_000001.11    +       65565   69037
XP_024307731.1  NC_000001.11    -       358067  399041
XP_024307730.1  NC_000001.11    -       358153  399041
NP_001005221.2  NC_000001.11    -       450740  450740
XP_011540840.1  NC_000001.11    -       586839  611112
```

## Configuring your Oxford Dot Plots

The plots produc
