<img src="https://github.com/conchoecia/odp/blob/main/figs/dotplot_fig_v2-01.jpg" width="250">

# odp - Oxford Dot Plots

## <a name="started"></a>Getting Started
```sh
git clone https://github.com/conchoecia/odp.git
# make a config.yaml file for your odp analysis
cp odp/example_configs/CONFIG_odp.yaml ./
# modify the config file to include your own data
vim CONFIG_odp.yaml
# rename the config so snakemake can find it
mv CONFIG_odp.yaml config.yaml
# run the pipeline
snakemake -r -p --snakefile odp/scripts/odp
# currently there is no man page, see https://github.com/conchoecia/odp/ for instructions
```

## Table of Contents

- [Getting Started](#started)
- [Users' Guide](#uguide)
  - [Installation](#install)
    - [Python Requirements](#python)
    - [Other Requirements](#otherreq)
  - [General usage](#general)
    - [Input File Requirements](#inputspec)
    - [`.chrom` file specifications](#chromspec)
    - [Help generating `.chrom` files](#chromhelp)
  - [Use cases](#cases)
    - [Make macrosynteny plots between two or more genomes](#macrosynuse)
    - [Find and characterize ancestral linkage groups](#ALGanalysis)
      - [ALGs part 1 - Ortholog finding in 3+ species](#nwayreciprocalbest)
      - [ALGs part 2 - Find significantly numerous groups of orthologs](#groupby)




## <a name="uguide"></a>Users' Guide

Odp is a protein-based synteny analysis software suite that is useful for comparing the evolution of chromosomes between two or more species. Use cases include (1) ploting synteny relationships between two genome assemblies, (2) inferring evolutionary relationships using chromosome synteny information, and (3) determining groups of ancestrally linked genes given a set of species' chromosome-scale genomes.

This software was visually modelled on the dotplots found in [Simakov, Oleg, et al. "Deeply conserved synteny resolves early events in vertebrate evolution." Nature ecology & evolution 4.6 (2020): 820-830.](https://www.nature.com/articles/s41559-020-1156-z), and was further expanded to determine the phylogenetic tree toplogy of animals in [Schultz, D.T., et al. (2022)](https://www.biorxiv.org/).

This software fills a niche in that it automates comparisons of chromosome-scale genomes, an increasingly important task as the genomes of more non-model organisms become available.

## <a name="install"></a>Installation

Odp and its dependencies are developed for a unix environment (linux, Mac OS X) running bash as the shell. You can download the software with this command:

```
git clone https://github.com/conchoecia/odp.git
```

### <a name="python"></a>Python Requirements

Your active python environment must be python 3. This software is implemented in [`snakemake`](https://snakemake.readthedocs.io/en/stable/). Specific python packages within the pipeline that must be available in your path are:

```
snakemake
matplotlib
pandas
numpy
```

If you have `conda` I recommend `conda install snakemake matplotlib pandas numpy` if you are not sure if you have the required packages.

### <a name="otherreq"></a>Other Requirements

Direct calls to these programs must also be available in your environment. Future versions of `odp` may bundle these software packages directly to avoid these requirements.

- [diamond](https://github.com/bbuchfink/diamond)
- awk
- [pauvre](https://github.com/conchoecia/pauvre) - I may remove this requirement in future releases 
- [bioawk](https://github.com/lh3/bioawk) - I may remove this requirement in future releases


For the aims above, this software works by:
1. For comparisons between two species, this program finds reciprocal-best protein matches using diamond blastp. The pipeline performs comparions between all *n* species in the config file. Compute time scales quadratically with increasing species *O(n*<sup>2</sup>*)*. The [Da and Db](https://www.nature.com/articles/s41559-020-1156-z) of each pairwise comparison is calculated to determine synteny block cutoffs in the cases of complex rearrangements. The signifiance of interactions between two or more genomes is calculated 
2. Calculating [Da and Db](https://www.nature.com/articles/s41559-020-1156-z) for each genome.
3. Plotting the genome assembly, reciprocal best protein hits, and Da/Db using matplotlib

## <a name="general"></a>General Usage

Odp requires, at minimum, the genome assembly sequence file, a sequence file of proteins found in the genome, and a file specifying the protein coordinates in the genome. The paths to these files for each sample is specified in a `.yaml` configuration file. 

A minimal working example of a config file that is set up to compare the genomes of _C. elegans_ and _H. sapiens_ looks like this:

```yaml
# this file is called config.yaml
xaxisspecies:
  Celegans:
    proteins: /path/to/proteins_in_Cel_genome.fasta
    prot_to_loc: /path/to/Cel_genome_annotation.chrom
    genome: /path/to/Cel_genome_assembly.fasta
  Homosapiens:
    proteins: /path/to/Human_prots.fasta
    prot_to_loc: /path/to/Human_annotation.chrom
    genome: /path/to/Human_genome_assembly.fasta
```

You can perform a comparison between these two genomes with:

```
snakemake --snakefile odp/scripts/odp
```

### <a name="inputspec"></a>Input file requirements

The file formats that are needed for the three files per genome are:
1. A genome assembly in [`.fasta` format](https://en.wikipedia.org/wiki/FASTA_format).
2. Protein sequences in [`.fasta` format](https://en.wikipedia.org/wiki/FASTA_format).
3. A file which details where the proteins are located in the genome, in [`.chrom` format](#chromspec)

### <a name="chromspec"></a>`.chrom` file specifications

The `.chrom` file format has 5 tab-delimited fields. Each line details the location of a protein on a scaffold. The fields are: 
  - `protein_header scaffold_header strand start stop`.

The requirements for each field are:
1. `protein_header` - the string here must match the header of a protein in the protein fasta.
2. `scaffold_header` - the string here must match the header of a sequence in the genome assembly fasta.
3. `strand` - must be `+` or `-`.
4. `start` - the numerically least position, in basepair coordinates, of the CDS of the protein. Like the start coords in a `bed` or `GFF` file, not necessarily the position of the start codon. Can often be found in a GFF3 or GTF.
5. `stop` - same as #4, but the stop position.

For example, the following `.chrom` file details four proteins that exist on two scaffolds. Two of the proteins are on the negative strand. The first protein, `BFGX8636T1`, has its start codon from the first position of scaffold 1, and the last codon ends at base 1246.

```
BFGX8636T1      sca1    +       1       1246
BFGX0001T1      sca1    -       2059    2719
BFGX0002T1      sca2    +       6491    12359
BFGX0003T1      sca2    -       12899   18848
```

### <a name="chromhelp"></a>Help generating `.chrom` files

A `.chrom` file can usually easily be generated from a genome annotation, such as a [`GFF3`](https://m.ensembl.org/info/website/upload/gff3.html) or [`GTF/GFF2`](https://uswest.ensembl.org/info/website/upload/gff.html) file. If you are working with NCBI GFFs, CDS entries have a predictable format that enables us to compile all of the information required for a chrom file: `NC_000001.11	BestRefSeq	CDS	65565	65573	.	+	0	ID=cds-NP_001005484.2;Parent=rna-NM_001005484.2;Dbxref=CCDS:CCDS30547.1,Ensembl:ENSP00000493376.2,GeneID:79501,Genbank:NP_001005484.2,HGNC:HGNC:14825;Name=NP_001005484.2;gbkey=CDS;gene=OR4F5;product=olfactory receptor 4F5;protein_id=NP_001005484.2;tag=MANE Select`. There are many CDS lines per gene, so a special parsing program is required.

The program bundled with `odp`, [scripts/NCBIgff2chrom.py](scripts/NCBIgff2chrom.py), parses gzipped/uncompressed NCBI GFFs and gets the full protein span in the genome. Running [scripts/NCBIgff2chrom.py](scripts/NCBIgff2chrom.py) on the [human GFF from NCBI](https://www.ncbi.nlm.nih.gov/genome/?term=human) with the command `python NCBIgff2chrom.py GCF_000001405.39_GRCh38.p13_genomic.gff.gz` results in a legal `.chrom` file with all of the proteins from the annotation. This file can be easily filtered later on.

```
NP_001005484.2  NC_000001.11    +       65565   69037
XP_024307731.1  NC_000001.11    -       358067  399041
XP_024307730.1  NC_000001.11    -       358153  399041
NP_001005221.2  NC_000001.11    -       450740  450740
XP_011540840.1  NC_000001.11    -       586839  611112
```

## <a name="cases"></a>Use cases

### <a name="macrosynuse"></a>Make macrosynteny plots between two or more genomes

Program: `odp/scripts/odp`
Input: `config.yaml` file with the species you wish to compare.
Output:
  - Tables of reciprocal best blastp hits between species.
  - PDF figures of Oxford dot plots comparing genome macrosynteny.
  - PDF figures of plots showing the degree of homology between chromosomes.
  - Oxford dot plots colored by predefined gene IDs, or by orthologs in other species

`config.yaml` format for running `odp/scripts/odp`:

```yaml
# this file is called config.yaml

ignore_autobreaks: True  # Can be True or False. True results in faster jobs
diamond_or_blastp: blastp # pick your protein search program. blastp is more sensitive.

prot_to_color:    # This is an optional parameter. odp_groupby_to_rbh generates a file 
  MyColorGroups1: /path/to/file_of_color_groups1.tsv        # that can be used here to
  ColorGroups2:   /path/to/another_color_group_file.tsv     # color proteins by groups
  
xaxisspecies:
  Celegans:
    proteins: /path/to/proteins_in_Cel_genome.fasta # required field
    prot_to_loc: /path/to/Cel_annot.chrom           # required field
    genome: /path/to/Cel_genome_assembly.fasta      # required field
    
    genus: "Caenorhabditis" # This is an optional field
    species: "elegans" # This is an optional field
    
    minscafsize: 1000000 # optional field. Sets minimum scaffold size to plot
    
    manual_breaks:    # optional field, tells the software to treat breaks
      - "I:50000"     #  as separate units for calculating the homology p-values
      - "IV:9000000"  #  with Fisher's exact test. Useful for plotting centromeres.
      - "II:99009"    #  Here, we tell the software that Cel chroms I, IV, II have breaks.
      
    chrom_to_color:                     # This is also an optional field, it tells odp
      "I:all": "#e11f26"                #  to plot different chromosomes' proteins with
      "II:0-50000": "#8d4b68"           #  different colors. This plots chromosome II
      "II:50000-9999999": "#3a7eb5"     #  with two colors and chromosome I with one.
      
    plotorder:    # This optional field tells the software to only plot the scaffolds
      - "I"       #  listed here, and to do it in this order. This is useful for plotting
      - "II"      #  comparisons between two species where you want a specific order for
      - "III"     #  both species.
    
  Homosapiens:
    proteins: /path/to/Human_prots.fasta
    prot_to_loc: /path/to/Human_annotation.chrom
    genome: /path/to/Human_genome_assembly.fasta
```

Run the pipeline with the command `snakemake -r -p --snakefile odp/scripts/odp`. The output files will be located in the folder `synteny_analysis/`. In this folder there are these folders:
  - `blastp_results/`
    - Files pertaining to the blastp searches between all the genomes. Also has files with the reciprocal best hits.
  - `db/`
    - Blast databases for searches.
  - `dvalue_table/`
    - Tables with the best reciprocal hits between species, genome coordinates.
  - `dvalue_table_breaks/`
    - The same files as in `dvalue_table/`, except that it contains *Da* and *Db* values for each pair.
  - `genome_coords/`
    - The genome coordinates (only used for plotting).
  - `odp_plotting_configs/`
    - Temporarily defunct, but generates `config.yaml` files that can be used for later steps in `odp`.
  - `plot_order/`
    - For each comparison, has the order in which scaffolds will be plotted.
  - `plots/`
    - `sample_similarity/`
      - Plots of distributions of ortholog length vs protein percent identity.
    - `significance/`
      - Plots containing several visual styles showing the significance of interactions between chromosomes, or sub-chromosomal regions defined by the user.
    - `synteny_colored_by*/`
      - Oxford dot plots of two species, color by the specified color from a list of protein ids and colors.
    - `synteny_colored_by_no_missing/`
      - Oxford dot plots of two species colored by the orthology and color identity of another species.
    - `synteny_uncolored`
      - Oxford dot plots. Every orthologous protein is colored blue.
  - `prot_to_color/`
    - A list of protein ids and how to color them.

### <a name="ALGanalysis"></a>Find and characterize ancestral linkage groups

Finding ancestral linkage groups of proteins for a group of species is a useful way to characterize what the genome at the ancestral node of that clade may have looked like, and to analyze how the genomes have evolved since that node.  See [Simakov et al. (2022)](https://www.science.org/doi/full/10.1126/sciadv.abi5884) for an example on how this concept was used to determine the ancestral number of chromosomal linkage groups in the common ancestor of sponges, cnidarians, and bilaterians.

The current implementation of this pipeline uses multiple steps to perform these analyses and determine the ALGs. For future versions of odp, we plan to implement this analysis into a single step.

#### <a name="nwayreciprocalbest"></a>ALGs part 1 - Ortholog finding in 3+ species

For this analysis, `blastp` or `diamond` analyses are performed against _n_ species that you specify. Orthologs are kept only when proteins in the _n_ speices are reciprocal best hits of each other. These are found by loading the `blast` results into a graph structure and finding [bidirectional complete graphs](https://en.wikipedia.org/wiki/Complete_graph) of `blastp` hits. This process is highly conservative, therefore as the number of genomes _n_ increases, the number of highly conserved orthologs decreases.

Program: `odp/scripts/odp_nway_rbh`
Input: `config.yaml`, same as for `odp/scripts/odp`, but with some modifications.
Output:
  - a file that contains the reciprocal best hits for the species included. This is called a `.rbh` file

`config.yaml` format for running `odp/scripts/odp`:

```yaml
# this file is called config.yaml

# the number of species you want to be included in each analysis
nways: 3
# How you want to identify the orthologs [diamond|blastp]
search_method: diamond
# What analyses you want to produce. Saves on some compute.
#  Must match headers of `xaxisspecies`. Order doesn't matter.
analyses:
  - ["Celegans", "Homosapiens", "Dmel"]
  - ["Celegans", "Homosapiens", "Mmus"]
  - ["Celegans", "Dmel", "Mmus"]
  # - ["Homosapiens", "Dmel", "Mmus"]   # You can comment out lines if you would like
  
xaxisspecies:
  Celegans:
    proteins: /path/to/proteins_in_Cel_genome.fasta
    prot_to_loc: /path/to/Cel_genome_annotation.chrom
    genome: /path/to/Cel_genome_assembly.fasta
  Homosapiens:
    proteins: /path/to/Human_prots.fasta
    prot_to_loc: /path/to/Human_annotation.chrom
    genome: /path/to/Human_genome_assembly.fasta
  Dmel:
    proteins: /path/to/drosophila_prots.fasta
    prot_to_loc: /path/to/drosophila_annotation.chrom
    genome: /path/to/drosophila_genome_assembly.fasta
  Mmus:
    proteins: /path/to/mouse_prots.fasta
    prot_to_loc: /path/to/mouse_annotation.chrom
    genome: /path/to/mouse_genome_assembly.fasta
```

The results of these analyses are found in `odp_nway_rbh/rbh/`. The reciprocal best hits file contains an unspecified number of columns, but always contains the columns:

* A unique identifier for each ortholog
* The protein ID for each species in that ortholog
* The scaffold ID on which that chromosome resides in that species
* The scaffold coordinates on which that chromosome resides in that species

Currently, the naming convention for these files is `[Sp.1]_[Sp.2]...[Sp.N]_reciprocal_best_hits.rbh`. This format is used in downstream steps to parse the headers found in the file.

#### <a name="groupby"></a>ALGs part 2 - Find significantly numerous groups of orthologs

The output of the previous program, the `.rbh` file, has one ortholog per line. In this step, we will group the orthologs together based on whether they exist on the same set of scaffolds in each species. For example, all of the orthologs that exist on:
  * _C. elegans_ chromosome I
  * _D. melanogaster_ chromosome 3
  * and human *chromosome 17*

will be one group. All of the orthologs that exist on a slightly different set of chromosomes will be another group, for example:  
  * _C. elegans_ chromosome I
  * _D. melanogaster_ chromosome 3
  * and *human chromosome 16*

The groups are saved in a tab-delimited file called a `.groupby` file. Each line is one group, and gene ids, scaffolds on which they reside, and genome coordinates are saved in python-type lists in single columns.

The number of groups found in this analysis, and the number of genes found in each group, depend on the degree of shared macrosynteny between the species used in the analysis. Distantly related species, or species with fast-evolving genomes, will have many groups, each with few genes. Closely related species or species with slowly-evolving genomes will have fewer groups, with more genes per group. Regardless of the relationships between the species, there will be a log decay of group sizes given the phenomenon of single genes translocating to other chromosomes.

For each group of orthologs *G*, we can estimate the false discovery rate *α* of finding a group of that size *|G|* given the specific genomes *s* in the analysis. We estimate this by taking the genomes *s*, and producing randomized versions of those genomes *s<sub>r</sub>* by shuffling the gene IDs in the `.chrom` file. This is performed millions of times, and we count the frequency of finding groups of size *|G|*. In other words, *α = 
