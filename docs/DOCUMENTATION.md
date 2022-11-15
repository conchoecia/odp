## Table of Contents

- [Introduction](#intro)
- [Specifications](#specs)
  - [LG databases](#lgdatabase)
    - [File 1 - the `.rbh` file](#LGrbhfile)
    - [File 2 - the `.hmm` file](#LGhmmfile)
- [Citations](#citations)


## <a name="intro"></a>Introduction

This goal of this documentation is to provide clear guidance on the use
and extension of this software package. Below, we provide specifications
on certain file formats.

## <a name="specs"></a>Specifications

Below are the specifications for a variety of file formats used in these analyses.

### <a name="lgdatabase"></a>LG databases

The LG databases provide a way to look at the relationships between species of
interest and conserved linkage groups of genes. For example, there are conserved
ancestral linkage groups of genes common to sponges, cnidarians, and bilaterians [1][sa2022].

The LG database consists of two files in a **one directory**.
  - File 1 is a special `.rbh` file that contains information about which orthologs are present in which linkage groups
  - File 2 is a `.hmm` with one HMM per ortholog in that is in the `.rbh` file

### <a name="LGrbhfile"></a>File 1 - the rbh file. (a very small example shown below).

Here, we see 9 orthologs (rows) in three linkage groups (`A1a`. `A2`. `B1`), and two species (`BFL`, `EMU`).

```
+──────────────+─────────────+──────────+───────────+─────────────────+──────────+──────────────+───────────+──────────+
| rbh          | gene_group  | color    | BFL_scaf  | BFL_gene        | BFL_pos  | EMU_gene     | EMU_scaf  | EMU_pos  |
+──────────────+─────────────+──────────+───────────+─────────────────+──────────+──────────────+───────────+──────────+
| BFL_EMU_134  | A1a         | #C23D51  | BFL1      | XP_035675733.1  | 13421192 | Em0019g254a  | EMU19     | 1407261  |
| BFL_EMU_135  | A1a         | #C23D51  | BFL1      | XP_035671256.1  | 16112138 | Em0019g644a  | EMU19     | 3783090  |
| BFL_EMU_136  | A1a         | #C23D51  | BFL1      | XP_035681411.1  | 27614506 | Em0019g708a  | EMU19     | 4077394  |
| BFL_EMU_137  | A1a         | #C23D51  | BFL1      | XP_035697917.1  | 31333283 | Em0019g141a  | EMU19     | 789406   |
| BFL_EMU_2177 | A2          | #8b4d68  | BFL1      | XP_035678153.1  | 23733754 | Em0006g1172a | EMU6      | 10781165 |
| BFL_EMU_2181 | A2          | #8b4d68  | BFL1      | XP_035672737.1  | 7868083  | Em0006g895a  | EMU6      | 8433533  |
| BFL_EMU_2186 | A2          | #8b4d68  | BFL1      | XP_035677008.1  | 4867558  | Em0006g1245a | EMU6      | 11271735 |
| BFL_EMU_2187 | A2          | #8b4d68  | BFL1      | XP_035686175.1  | 4661712  | Em0006g1244a | EMU6      | 11266662 |
| BFL_EMU_57   | B1          | #fffb58  | BFL10     | XP_035688512.1  | 7539304  | Em0004g1007a | EMU4      | 8057909  |
+──────────────+─────────────+──────────+───────────+─────────────────+──────────+──────────────+───────────+──────────+
```

The `.rbh` file must have the columns: `["rbh", "gene_group", "color"]`. Each row represents one ortholog. The first row is the headers.

- Columns:
  - `rbh` - is an arbitrary name of the ortholog. Must be unique within this file. Doesn't need to have the species' names 
  - `gene_group` - is the group to which this ortholog belongs. These letters are from [1][sa2022].
  - `color` - is the color of that dot. We don't enforce that all entries of a group must share the same color, so proceed with caution if you wish to be consistent within a group.

- Additional Columns:
  - For each species that is present in the `.rbh` file there must be three columns. The prefix of the three additional columns is the species name. The species name must not contain a `_` character.
    - `Species_gene` - the protein fasta header of this this species' gene in the ortholog
    - `Species_scaf` - the scaffold on which the protein exists in this species' genome
    - `Species_pos` - the position of the protein on its scaffold.

### <a name="LGhmmfile"></a>File 2 - the hmm file.

The directory containing the `.rbh` file must also contain a protein `.hmm` file in the HMMER format from the Eddy lab [2][hmmer].
Within the `.hmm` file, there must be one HMM from the output of `hmmbuild` for each ortholog in the `rbh` file.
The name of the HMM must be the same name as in the `rbh` column of the `.rbh` file. For example, the entry
for the ortholog `BFL_EMU_134` would start like this in the `.hmm` file:

```
HMMER3/f [3.3.2 | Nov 2020]
NAME BFL_EMU_134
.
.
.
```

One must build the alignments by aligning the proteins within each ortholog. We use `mafft` in several `odp` scripts, however
other alignment methods are fine. With the alignment, each `.hmm` must be built with `hmmbuild`.

## <a name="citations"></a>Citations
  1. [Simakov, O., Bredeson, J., Berkoff, K., Marletaz, F., Mitros, T., Schultz, D.T., O’Connell, B.L., Dear, P., Martinez, D.E., Steele, R.E. and Green, R.E., 2022. Deeply conserved synteny and the evolution of metazoan chromosomes. Science advances, 8(5), p.eabi5884.][sa2022]
  2. [Eddy, S.R., 1995, July. Multiple alignment using hidden Markov models. In Ismb (Vol. 3, pp. 114-120).][hmmer]

[sa2022]: https://www.science.org/doi/10.1126/sciadv.abi5884
[hmmer]: http://hmmer.org/
