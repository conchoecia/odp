<img src="https://github.com/conchoecia/odp/blob/main/docs/dotplot_fig_v2-01.jpg" width="250">

# odp - Oxford dot plots

## <a name="started"></a>Getting Started
```sh
#install
git clone https://github.com/conchoecia/odp.git
cd odp && make
# make a config.yaml file for your odp analysis
cp odp/example_configs/CONFIG_odp.yaml ./config.yaml
# modify the config file to include your own data
vim config.yaml
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
      - [ALGs part 3 - Filter groups of orthologs](#groupbyfilter)  
      - [ALGs part 4 - Annotate groups of orthologs](#groupbyannotate)
      - [ALGs part 5 - Merge groupby files](#groupbymerge) 
      - [ALGs part 6 - Find orthologs in more species](#groupbytohmm)
      - [ALGs part 7 - Plot mixing of select linkage groups](#plotmixing)
    - [Determine which clade is sister](#4speciesphylogeny)

## <a name="uguide"></a>Users' Guide

Odp is a protein-based synteny analysis software suite that is useful for
comparing the evolution of chromosomes between two or more species. Use cases
include (1) ploting synteny relationships between two genome assemblies, (2)
inferring evolutionary relationships using chromosome synteny information, and
(3) determining groups of ancestrally linked genes given a set of species'
chromosome-scale genomes.

This software was visually modelled on the dotplots found in [Simakov, Oleg, et al. "Deeply conserved synteny resolves early events in vertebrate evolution." Nature ecology & evolution 4.6 (2020):820-830.](https://www.nature.com/articles/s41559-020-1156-z), and was further
expanded to determine the phylogenetic tree toplogy of animals in
[Schultz, D.T., et al. (2023)](https://www.biorxiv.org/).

This software fills a niche in that it automates comparisons of chromosome-scale
genomes, an increasingly important task as the genomes of more non-model
organisms become available.

## <a name="install"></a>Installation

Odp and its dependencies are developed for a unix environment (linux, Mac OS X)
running bash as the shell. You can download the software with this command:

```
git clone https://github.com/conchoecia/odp.git
cd odp && make
```

### <a name="python"></a>Python Requirements

Your active python environment must be python 3. This software is implemented in
[`snakemake`](https://snakemake.readthedocs.io/en/stable/). Specific python
packages within the pipeline that must be available in your python installation
are:

```
snakemake
matplotlib
networkx
scipy
pandas
numpy
seaborn
```

If you have `conda` I recommend `conda install snakemake matplotlib pandas numpy seaborn`
if you are not sure if you have the required packages.

### <a name="otherreq"></a>Other Requirements

Direct calls to these programs must also be available in your environment.
Future versions of `odp` may bundle these software packages directly to avoid
these requirements.

- [diamond](https://github.com/bbuchfink/diamond)
- [blastp](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
- awk

For the aims above, this software works by:
1. For comparisons between two species, this program finds reciprocal-best protein matches using diamond blastp. The pipeline performs comparions between all *n* species in the config file. Compute time scales quadratically with increasing species *O(n*<sup>2</sup>*)*. The [Da and Db](https://www.nature.com/articles/s41559-020-1156-z) of each pairwise comparison is calculated to determine synteny block cutoffs in the cases of complex rearrangements. The signifiance of interactions between two or more genomes is calculated 
2. Calculating [Da and Db](https://www.nature.com/articles/s41559-020-1156-z) for each genome.
3. Finding which orthologs correspond to previously identified ancestral linkage groups (ALGs), such as the Bilateria-Cnidaria-Sponge ALGs from [Simakov et al. (2020).](https://www.nature.com/articles/s41559-020-1156-z).
3. Plotting the genome assembly, reciprocal best protein hits, the ALGs, and Da/Db using matplotlib.

## <a name="general"></a>General Usage

Odp requires, at minimum, the genome assembly sequence file, a sequence file of proteins found in the genome, and a file specifying the protein coordinates in the genome. The paths to these files for each sample is specified in a `.yaml` configuration file. 

A minimal working example of a config file that is set up to compare the genomes of _C. elegans_ and _H. sapiens_ looks like this:

```yaml
# this file is called config.yaml
ignore_autobreaks: True       # Skip steps to find breaks in synteny blocks
diamond_or_blastp: "diamond"  # "diamond" or "blastp"
plot_LGs: True                # Plot the ALGs based on the installed databases
plot_sp_sp: True              # Plot the synteny between two species, if False just generates .rbh files

species:
  Celegans:
    proteins: /path/to/proteins_in_Cel_genome.fasta
    chrom: /path/to/Cel_genome_annotation.chrom
    genome: /path/to/Cel_genome_assembly.fasta
  Homosapiens:
    proteins: /path/to/Human_prots.fasta
    chrom: /path/to/Human_annotation.chrom
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
3. A file which details where the proteins are located in the genome, in [`.chrom` format](#chromspec). Using a `.gff` or `.gtf` file is currently not supported, but support is planned.

### <a name="chromspec"></a>`.chrom` file specifications

The `.chrom` file format has 5 tab-delimited fields. Each line details the location of a protein on a scaffold. The fields are:
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
  
species:
  Celegans:
    proteins: /path/to/proteins_in_Cel_genome.fasta # required field
    chrom:    /path/to/Cel_annot.chrom              # required field
    genome:   /path/to/Cel_genome_assembly.fasta    # required field 

    genus: "Caenorhabditis" # This is an optional field
    species: "elegans" # This is an optional field 

    minscafsize: 1000000 # optional field. Sets minimum scaffold size to plot. 

    manual_breaks:    # optional field, tells the software to treat breaks
      - "I:50000"     #  as separate units for calculating the homology p-values
      - "IV:9000000"  #  with Fisher's exact test. Useful for plotting centromeres.
      - "II:99009"    #  Here, we tell the software that Cel chroms I, IV, II have breaks.      

    plotorder:    # This optional field tells the software to only plot the scaffolds
      - "I"       #  listed here, and to do it in this order. This is useful for plotting
      - "II"      #  comparisons between two species where you want a specific order for
      - "III"     #  both species.
    
  Homosapiens:
    proteins: /path/to/Human_prots.fasta
    chrom:    /path/to/Human_annotation.chrom
    genome:   /path/to/Human_genome_assembly.fasta
```

Run the pipeline with the command `snakemake -r -p --snakefile odp/scripts/odp`. The output files will be located in the folder `synteny_analysis/`. In this folder there are these folders:
  - `db`
    - blastp and and diamond databases for searches.
  - `step0-blastp_results`
    - blastp/diamond searches between all the genomes. Also has files with the reciprocal best hits.
  - `step0-chromsize`
    - chromosome information used later for plotting.
  - `step1-rbh`
    - reciprocal best hits files for the analyses.
  - `step2-figures`
    - `ALG-species_plots`
      - If there are ALGs installed in the `LG_db`, then the ALG-species plots will appear here.
    - `synteny_coloredby_*`
      - If there are ALGs installed in the `LG_db`, then the species-species synteny plots will appear here.
    - `synteny_nocolor`
      - Two-species synteny plots appear here regardless of what is in `LG_db`.

### <a name="ALGanalysis"></a>Find and characterize ancestral linkage groups

Finding ancestral linkage groups of proteins for a group of species is a useful
way to characterize what the genome at the ancestral node of that clade may have
looked like, and to analyze how the genomes have evolved since that node.  See
[Simakov et al. (2022)](https://www.science.org/doi/full/10.1126/sciadv.abi5884)
for an example on how this concept was used to determine the ancestral number of
chromosomal linkage groups in the common ancestor of sponges, cnidarians, and
bilaterians.

The current implementation of this pipeline uses multiple steps to perform these
analyses and determine the ALGs. For future versions of odp, we plan to
implement this analysis into a single step.

#### <a name="nwayreciprocalbest"></a>ALGs part 1 - Ortholog finding in 3+ species

For this analysis, `blastp` or `diamond` analyses are performed against _n_
species that you specify. Orthologs are kept only when proteins in the _n_
speices are reciprocal best hits of each other. These are found by loading the
`blast` results into a graph structure and finding
[bidirectional complete graphs](https://en.wikipedia.org/wiki/Complete_graph) of `blastp`
hits. This process is highly conservative, therefore as the number of genomes _n_
increases, the number of highly conserved orthologs decreases.

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
  
species:
  Celegans:
    proteins: /path/to/proteins_in_Cel_genome.fasta
    chrom: /path/to/Cel_genome_annotation.chrom
    genome: /path/to/Cel_genome_assembly.fasta
  Homosapiens:
    proteins: /path/to/Human_prots.fasta
    chrom: /path/to/Human_annotation.chrom
    genome: /path/to/Human_genome_assembly.fasta
  Dmel:
    proteins: /path/to/drosophila_prots.fasta
    chrom: /path/to/drosophila_annotation.chrom
    genome: /path/to/drosophila_genome_assembly.fasta
  Mmus:
    proteins: /path/to/mouse_prots.fasta
    chrom: /path/to/mouse_annotation.chrom
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

The groups are saved in a tab-delimited file called a `.groupby` file. Each line
is one group, and gene ids, scaffolds on which they reside, and genome
coordinates are saved in python-type lists in single columns.

The number of groups found in this analysis, and the number of genes found in
each group, depend on the degree of shared macrosynteny between the species used
in the analysis. Distantly related species, or species with fast-evolving
genomes, will have many groups, each with few genes. Closely related species or
species with slowly-evolving genomes will have fewer groups, with more genes per
group. Regardless of the relationships between the species, there will be a log
decay of group sizes given the phenomenon of single genes translocating to other
chromosomes.

For each group of orthologs *G*, we can estimate the false discovery rate *α* of
finding a group with *i* or fewer genes given the real genomes in these
comparisons. We estimate this false discovery rate by producing randomized
versions of the genomes by shuffling the gene IDs in the `.chrom` file,
measuring whether a group of *i* or fewer genes was present, then repeating this
measurement hundreds of millions of times.

#### <a name=“groupbyfilter”></a>ALGs part 3 - Filter groups of orthologs

The groups of reciprocal best hits, as well as the newly-calculated false
discovery rates, are saved in the resulting `.groupby` file. This can be
manually or programmatically filtered to only keep groups with certain
properties, or groups with a significantly low false discovery rate.

This is performed automatically by `odp_nway_rbh`, but can be performed with the
script `odp_groupby_filter` by specifying the `.groupby` file, and by specifying
the acceptable false discovery rate cutoff.

This process can also be performed in a table editor, such as the spreadsheets
on Google Drive, Apple Sheets, or Microsoft Excel.

After removing the rows that have a less-than-significant false discovery rate,
continue on to the next step to annotate the groups of orthologs.

#### <a name=“groupbyannotate”></a>ALGs part 4 - Annotate groups of orthologs

At this stage the resulting rows are groups of orthologous genes that are
present on the same set of chromosomes in the species under consideration, and
have been since the common ancestor of these species. In other words, these are
ancestral linkage groups (ALGs) for this clade.

It is useful at this stage to assign names to each of the rows in the `group` column of the `.groupby` file. There is some precedence for these naming conventions, see [Simakov et al. (2020)](https://www.nature.com/articles/s41559-020-1156-z) and [Simakov et al. (2022)](https://www.science.org/doi/full/10.1126/sciadv.abi5884). So, if your analysis includes animal genomes then it may be helpful to include some of the species from these publications.

It is not necessary that each row has its own unique group ID. However, doing so will help plot mixing in downstream analyses.

#### <a name=“groupbymerge”></a>ALGs part 5 - Merge `.groupby`/`.rbh` files

In this section let’s consider a few species to compare:
  - Unicellular (+colonial multicellular) outgroups of animals:
    - the icthyosporean _Creolimax fragrantissima_ (CFR)
    - the filasterean amoeba _Capsaspora owczarzaki_ (COW)
    - the choanoflagellate _Monosiga brevicollis_ (MBR)
  - The animals:
    - the ctenophore _Hormiphora californensis_ (HCO)
    - the sponge _Ephydatia muelleri_ (EMU)
    - the jellyfish _Rhopilema esculentum_ (RES)

It is desirable to merge the `.groupby` files from the searches of multiple species if the evolutionary distance between the outgroup and the other species is extreme. For example, the degree of synteny between animals and their unicellular Holozoan relatives is very little, and merging multiple searches enables the discovery of more ancestrally linked genes.

We can perform the ALG-finding steps 1-4 described above for the following three analyses:
  - `CFR-HCA-EMU-RES`
  - `COW-HCA-EMU-RES`
  - `MBR-HCA-EMU-RES`

Each row in the `.groupby` files for these analyses will contain one gene per species in the analysis. It is possible that many of the orthologs will also contain proteins in two or more unicellular outgroups, so we now run `odp_groupby_to_rbh` to unwrap each `.groupby` file to a `.rbh` file, then `odp_rbh_merge` to join the `.rbh` files on the species `HCA`, `EMU`, and `RES`.

Each ortholog (row) in the resulting `.rbh` file will have a gene for each animal species (`HCA`, `EMU`, `RES`), and will contain a gene in between one and three of the unicellular species (`CFR`, `COW`, `MBR`).

The notation we use to refer to an `.rbh` file created by merging other `.rbh` files uses parentheses to note the species that may have missing data, and unmodified text to note the species that will always have a gene for each ortholog. The analysis discussed above is notated as `(CFR-COW-MBR)-HCA-EMU-RES`.

#### <a name=“groupbytohmm”></a>ALGs part 6 - Find orthologs in more species

Steps 1-4 of finding ALGs relies on using only a few species (perhaps 3-5) to avoid loss of orthologs due to the stringent ortholog selection process. [Step 5 - Merge `.groupby`/`.rbh` files, discussed above,](#groupbymerge) enables the inclusion of more genes by allowing for missing data in select groups. Then, by constructing hidden Markov models of the orthologs, we can search for orthologs in more species.

The script `odp_rbh_to_hmm` reads in a `.rbh` file and constructs one HMM model per ortholog (row). The models are then searched against the proteins of every additional species that is included in the `config.yaml` file. The best protein for each HMM is selected, and only proteins with a significant match are kept. Missing data are permissible in this step, so it is not guaranteed that every ortholog will have an identifiable protein in every species added in this step.

The output of this pipeline is another `.rbh` file, now with the proteins of the additional species identified with the HMM.

#### <a name=“plotmixing”></a>ALGs part 7 - Plot mixing of select linkage groups

If you have followed the above steps, you now have a `.rbh` file with orthologs that have been annotated by group, and includes many additional species thanks to the merging and HMM search steps.

If you are using `odp` to look for phylogenetically diagnostic fusion-then-mixing events, then it is useful to plot linkage groups to visualize the extent of mixing of those groups. The script `odp_rbh_plot_mixing` does that. The output of this script are PNG and PDFs of the orthologs in those two groups plotted in the chromosome coordinates for each species. This script also estimates the degree of intermixing of two groups of genes on the chromosomes on which they coexist.

### <a name=“4speciesphylogeny”></a>Determine which clade is sister

The module `odp_genome_rearrangegment_simulation` was developed to help answer the question of whether ctenophores or sponges are the sister clade of all other animals. This script requires one species that is the known outgroup, and one species nested in the phylogeny with a known relative position to all other species. In our study, we performed analyses in which the filasterean amoeba _Capsaspora owczarzaki_ or the choanoflagellate _Salpingoeca rosetta_ were the outgroup species. The species with the known phylogenetic position was the fire jellyfish _Rhopilema esculentum_. The program uses these genomes to polarize the relationships between the two genomes in an unresolved polytomy, in this case the ctenophore _Hormiphora californensis_ and the sponge _Ephydatia muelleri_.

The program `odp_genome_rearrangement_simulation` does the following:
1. Finds linkage groups that simultaneously satisfy all of these requirements:
  - The linkage groups are on separate chromosomes in the outgroup species
  - The linkage groups are on separate chromosomes in one of the unplaced species
  - The linkage groups are fused and mixed on single chromosomes in the other unplaced species
  - The linkage groups are fused and mixed on single chromosomes in the species with a known phylogenetic positions
2. Quantifies the number of fusion events identified in step 1, and calculates the number of genes per linkage group, and the number of linkage groups participating in those events.
3. For each species in the analysis, tests whether the phylogentically informative fusion-then-mixing events seen in the real genomes are due to randomness, or true biological signal. This is done by:
  - One simulation shuffles the protein IDs in the genome of one species. Step 1 above is run on this shuffled genome, and the three other observed genomes. In other workds, the number of phylogenetically informative fusion-then-mixing events are identified given the new shuffled genome.
  - Step 2 above is run on the events found above.
  - These steps are performed millions of times to estimate how many times we see a genome configuration that has at least as many genes, linkage groups, and fusion events participating in phylogentically diagnostic fusion-then-mixing events.

The output of this program is histograms showing the different measured parameters from the simulations (grey bars), plotted with the parameters observed from the real genomes (red vertical bars). These plots show whether the sister clade hypotheses seen in the real data can be explained by a highly rearranged state in any of the genomes in the analysis.

## <a name=“citation”></a>Citation

There is currently no citation for this work, so please cite the repository. This section will be updated when the manuscript has been published or uploaded to bioRxiv.
