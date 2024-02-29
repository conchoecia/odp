"""
Program  : GenDB_build_db_unannotated.snakefile
Language : snakemake
Date     : 2023-10-01
Author   : Darrin T. Schultz
Email    : darrin.schultz@univie.ac.at
Github   : https://github.com/conchoecia/odp
License  : GNU GENERAL PUBLIC LICENSE, Version 3, 29 June 2007. See the LICENSE file.

Description:
  This program:
    - takes the list of unannotated genomes,
    - annotates the genomes with common animal linkage groups,
    - and prepares a database of those linkage groups for ODP.
  The program requires that either a directory of the annotated and unannotated genome lists be provided.
  Otherwise the user has to specify specific paths to those tsv files.

Usage instructions:
  - There are currently no usage instructions. This is a work in progress.
"""

# This block imports fasta-parser as fasta
import os
import sys
snakefile_path = os.path.dirname(os.path.realpath(workflow.snakefile))
dependencies_path = os.path.join(snakefile_path, "../dependencies/fasta-parser")
sys.path.insert(1, dependencies_path)
import fasta

import pandas as pd
from datetime import datetime
import GenDB
import yaml

# figure out where bin is because we need to use some outside tools
snakefile_path = os.path.dirname(os.path.realpath(workflow.snakefile))
bin_path = os.path.join(snakefile_path, "../bin")

#if "API_key" not in locals():
#    API_key = ""

configfile: "config.yaml"
config["tool"] = "odp_ncbi_genome_db"
# Do some logic to see if the user has procided enough information for us to analyse the genomes
if ("directory" not in config) and ("accession_tsvs" not in config):
    raise IOError("You must provide either a directory of the annotated and unannotated genome lists, or a list of the paths to those tsv files. Read the config file.")

config["tempdir"] = "/tmp"
# check that the tempdir exists
if "tempdir" not in config:
    raise IOError("You must provide a temporary directory to store temporary files. Read the config file.")
# strip all trailing slashes from the tempdir
config["temp"] = config["tempdir"].rstrip("/").rstrip("\\")
if not os.path.isdir(config["tempdir"]):
    raise IOError("The temporary directory you provided does not exist. {}".format(config["tempdir"]))

config = GenDB.opening_logic_GenDB_build_db(config, chr_scale = True, annotated = False)

# One key feature of this script is that we will map proteins from
#  LG databases to annotate those genomes with the LG identities.
#  Before we accept the user's input, we need to parse the supplied
#  directories to make sure that they are valid.
LG_to_db_directory_dict = {}
LG_to_rbh_dfs           = {}
LG_outfiles             = []
if len(LG_to_db_directory_dict) == 0: # only do this once
    if "LG_db_directories" not in config:
        raise ValueError("You have not specified the directories of LG databases. Please add the key 'LG_db_directories' to your config file.")
    for thisdirectory in config["LG_db_directories"]:
        # The directory must exist
        if not os.path.isdir(thisdirectory):
            raise IOError("The directory of LG databases you provided does not exist. {}".format(thisdirectory))
        # There must be a directory named aligned
        if not os.path.isdir(os.path.join(thisdirectory, "aligned")):
            raise IOError("The directory of LG databases you provided does not contain a directory named 'aligned'. {}".format(thisdirectory))
        # The directory name is the name of the LG, just get the last part of the path
        LG_name = os.path.basename(os.path.normpath(thisdirectory))
        # The directory must contain a .rbh file
        if not any([x.endswith(".rbh") for x in os.listdir(thisdirectory)]):
            raise IOError("The directory of LG databases you provided does not contain a .rbh file. {}".format(thisdirectory))
        # Now we read in the .rbh file as a pandas df to get the LG names from the `.rbh` column.
        rbhfilepath = [x for x in os.listdir(thisdirectory) if x.endswith(".rbh")][0]
        df = pd.read_csv(os.path.join(thisdirectory, rbhfilepath), sep="\t")
        # For each rbh entry in the dataframe, there should be a fasta file in the aligned directory
        fasta_list = [x for x in os.listdir(os.path.join(thisdirectory, "aligned")) if x.endswith(".fasta")]
        for rbh in df["rbh"].unique():
            if rbh + ".fasta" not in fasta_list:
                raise IOError("The directory of LG databases you provided does not contain a fasta file for the rbh entry {}. {}".format(rbh, thisdirectory))
        # This appears to be a legitimate LG database directory. Add it to the dictionary.
        # Also add a pandas dataframe to the dictionary.
        LG_to_db_directory_dict[LG_name] = thisdirectory
        LG_to_rbh_dfs[LG_name] = df
        outfile =config["tool"] + "/input/LG_proteins/{}.fasta".format(LG_name)
        LG_outfiles.append(outfile)

resource_scopes:
    mem_mb="local",
    time="local",
    download_slots="global"

wildcard_constraints:
    taxid="[0-9]+",

# first we must load in all of the files. Only do it once
rule all:
    input:
        expand(config["tool"] + "/output/source_data/unannotated_genomes/{assemAnn}/{assemAnn}.chr.fasta.gz", assemAnn=config["assemAnn"]),
        #LG_outfiles,
        expand(config["tool"] + "/output/source_data/unannotated_genomes/{assemAnn}/{assemAnn}_annotated_with_{LG_name}.chrom.gz",
               assemAnn=config["assemAnn"], LG_name=LG_to_db_directory_dict.keys()),
        expand(config["tool"] + "/output/source_data/unannotated_genomes/{assemAnn}/{assemAnn}_annotated_with_{LG_name}.pep.gz",
               assemAnn=config["assemAnn"], LG_name=LG_to_db_directory_dict.keys()),
        expand("NCBI_odp_db.unannotated.{LG_name}.yaml",
               LG_name=LG_to_db_directory_dict.keys()),
        expand("NCBI_odp_sp_list.unannotated.{LG_name}.txt",
               LG_name=LG_to_db_directory_dict.keys()),

def dlChrs_get_mem_mb(wildcards, attempt):
    """
    The amount of RAM needed for the script depends on the size of the input genome.
    """
    attemptdict = {1: 4000,
                   2: 16000,
                   3: 64000
                  }
    return attemptdict[attempt]

rule dlChrs:
    """
    We have selected the unannotated genomes to download.
    For these genomes we need to find a way to annotate them.

    To specifically download the chromosome-scale scaffolds, there is a series of commands with the NCBI datasets tool.
      I found this set of instructions after opening a ticket on NCBI's github page: https://github.com/ncbi/datasets/issues/298

    Notes:
      - 12-31-2023 - Using the command line had too many edge cases that didn't work, so I resorted to using python to do a more careful job.
        This verifies that the files are downloaded and unzipped correctly, and contain all of the expected sequences.
        Therefore, we do not need additional verification steps for the assembly file.
    """
    input:
        datasets = os.path.join(bin_path, "datasets")
    output:
       fasta   = temp(config["tool"] + "/output/source_data/unannotated_genomes/{assemAnn}/{assemAnn}.chr.fasta"),
       allscaf = config["tool"] + "/output/source_data/unannotated_genomes/{assemAnn}/{assemAnn}.scaffold_df.all.tsv",
       chrscaf = config["tool"] + "/output/source_data/unannotated_genomes/{assemAnn}/{assemAnn}.scaffold_df.chr.tsv"
    retries: 3
    params:
        outdir   = config["tool"] + "/output/source_data/unannotated_genomes/{assemAnn}/",
    threads: 1
    group: "dlgz"
    resources:
        mem_mb = dlChrs_get_mem_mb, # 1 GB of RAM
        time   = 20,  # 20 minutes.
        download_slots = 1
    run:
        result = GenDB.download_unzip_genome(wildcards.assemAnn, params.outdir,
                                             input.datasets, chrscale = True)
        if result != 0:
            raise ValueError("The download of the genome {} failed.".format(wildcards.assemAnn))

rule gzip_fasta_file:
    """
    In this step zip the fasta file to conserve space.
    """
    input:
        genome = config["tool"] + "/output/source_data/unannotated_genomes/{assemAnn}/{assemAnn}.chr.fasta",
    output:
        genome = config["tool"] + "/output/source_data/unannotated_genomes/{assemAnn}/{assemAnn}.chr.fasta.gz"
    threads: 1
    group: "dlgz"
    resources:
        mem_mb = 1000, # 1 GB of RAM
        time   = lambda wildcards: GenDB.gzip_get_time(config["assemAnn_to_scaflen"][wildcards.assemAnn])
    shell:
        """
        echo "Gzipping the fasta file."
        gzip < {input.genome} > {output.genome}
        """

rule generate_LG_fasta_sequence:
    """
    Currently, in the LG database, the sequences are only available as alignments.
    Here, we just concatenate all of the alignments, then strip the gaps.
    """
    input:
        LG_dir = lambda wildcards: LG_to_db_directory_dict[wildcards.LG_name]
    output:
        fasta = config["tool"] + "/input/LG_proteins/{LG_name}.fasta"
    resources:
        mem_mb = 1000, # 1 GB of RAM
        time   = 5     # 5 minutes
    threads: 1
    run:
        with open(output.fasta, "w") as o:
            for fastafile in os.listdir(input.LG_dir + "/aligned"):
                for record in fasta.parse(input.LG_dir + "/aligned/" + fastafile):
                    o.write(">{}\n{}\n".format(record.id, record.seq.replace("-", "")))

def miniprot_get_mem_mb(wildcards, attempt):
    """
    The amount of RAM needed for miniprot is highly variable.
    """
    attemptdict = {1: 16000,
                   2: 32000,
                   3: 64000,
                   4: 128000,
                   5: 256000,
                   6: 512000,
                   7: 1024000,
                   8: 1536000}
    return attemptdict[attempt]

rule miniprot:
    """
    This handles all of the miniprot steps. Both the indexing and the mapping.
    Doing it this way prevents keeping a ton of temporary data on the hard drive.
    """
    input:
        pep    = config["tool"] + "/input/LG_proteins/{LG_name}.fasta",
        genome = config["tool"] + "/output/source_data/unannotated_genomes/{assemAnn}/{assemAnn}.chr.fasta.gz",
    output:
        paf  = config["tool"] + "/output/mapped_reads/{assemAnn}/{LG_name}_to_{assemAnn}.paf"
    threads: 8
    retries: 8
    params:
        mpi_suffix = "{LG_name}_{assemAnn}.filt.fasta.gz.mpi" # this is the temporary index file
    resources:
        tmpdir = config["tempdir"], # the place where the temporary index file will be stored
        mem_mb = miniprot_get_mem_mb, # The RAM usage can blow up during indexing. Often > 10GB. 6Gbp genomes need more than 20GB of RAM.
        time   = 60 # 20 minutes
    shell:
        """
        INDEXFILE=$TMPDIR/{params.mpi_suffix}
        # don't index if the file already exists
        if [ ! -f ${{INDEXFILE}} ]; then
            # INDEXING step
            echo "Indexing the genome: Genome {wildcards.assemAnn} for LG {wildcards.LG_name}"
            miniprot -t {threads} -d ${{INDEXFILE}} {input.genome}
        fi

        # MAPPING step
        echo ""
        echo "Mapping the proteins to the genome. {wildcards.LG_name} to {wildcards.assemAnn}"
        miniprot -t {threads} ${{INDEXFILE}} {input.pep} > {output.paf}

        # REMOVE the index
        echo ""
        echo "Removing the index file: ${{INDEXFILE}}"
        rm ${{INDEXFILE}}
        """

rule filter_paf_for_longer_scaffold:
    """
    Many times, a single protein will map equally well to multiple scaffolds.
    In this case, we just pick the longest of the scaffolds.

    To do this, we load of a file of all the chrom sizes first to do the comparisons.
    """
    input:
        paf = config["tool"] + "/output/mapped_reads/{assemAnn}/{LG_name}_to_{assemAnn}.paf",
    output:
        paf = config["tool"] + "/output/mapped_reads/{assemAnn}/{LG_name}_to_{assemAnn}.filt.paf"
    threads: 1
    resources:
        mem_mb = 100, # I can't forsee using a GB of RAM, but easy to request.
        time   = 1    # Just 10 minutes
    run:
        paf_colnames = ["query",  "qlen", "qstart", "qend", "strand",
                         "target", "tlen", "tstart", "tend", "matches",
                         "alnlen", "mapq", "AS", "ms"]

        entries = []
        with open(input.paf, "r") as f:
            for line in f:
                # skip the header lines
                if line.startswith("#"):
                    continue
                # split the line
                line = line.strip().split("\t")[:len(paf_colnames)]
                line[-2] = line[-2].replace("AS:i:", "")
                line[-1] = line[-1].replace("ms:i:", "")
                for colindex in [1,2,3,6,7,8,9,10,11,12,13]:
                    line[colindex] = int(line[colindex])
                entries.append(line)
        # load in the paf file, only read in the first 14 columns
        pafdf = pd.DataFrame(entries, columns=paf_colnames)
        # Sort by query, then ms, then tlen. This will put the best hits at the top.
        # Then, drop duplicates of the query column. This will keep only the best hit.
        pafdf = pafdf.sort_values(by=["query", "ms", "tlen"], ascending=False
                 ).drop_duplicates(subset=["query"])
        # save the filtered paf file to the outout, but don't use the index or header
        pafdf.to_csv(output.paf, sep="\t", index=False, header=False)

rule paf_to_chrom_and_pep:
    """
    This script converts the filtered paf file to the chrom file.
    There are some decisions to make here about which proteins are best

    The chrom format is: protein_name, scaf_name, strand, scaf_start, scaf_end
    BFGX8636T1      sca1    +       1       1246
    BFGX0001T1      sca1    -       2059    2719
    BFGX0002T1      sca2    +       6491    12359
    BFGX0003T1      sca2    -       12899   18848

    Also outputs a protein fasta file of the entries.
    """
    input:
        paf = config["tool"] + "/output/mapped_reads/{assemAnn}/{LG_name}_to_{assemAnn}.filt.paf",
        pep  = config["tool"] + "/input/LG_proteins/{LG_name}.fasta"
    output:
        chrom = temp(config["tool"] + "/output/source_data/unannotated_genomes/{assemAnn}/{assemAnn}_annotated_with_{LG_name}.chrom"),
        pep   = temp(config["tool"] + "/output/source_data/unannotated_genomes/{assemAnn}/{assemAnn}_annotated_with_{LG_name}.pep")
    threads: 1
    resources:
        mem_mb = 500, # I can't forsee using a GB of RAM, but easy to request.
        time   = 5    # Just 5 minutes
    run:
        paf_colnames = ["query",  "qlen", "qstart", "qend", "strand",
                         "target", "tlen", "tstart", "tend", "matches",
                         "alnlen", "mapq", "AS", "ms"]

        # read in the paf file as a pandas dataframe
        df = pd.read_csv(input.paf, sep="\t", header=None, names=paf_colnames)
        # Split the query column into two parts. The first part is up to the last
        # _ character, the second part is everything else. Save the first part as
        # the rbh column.
        df["rbh"]     = df["query"].str.rsplit("_", n=1, expand=True)[0]
        # the second part is the species name
        df["species"] = df["query"].str.rsplit("_", n=1, expand=True)[1]
        # sort based on the rbh column, then the ms column, then the tlen column
        df = df.sort_values(by=["rbh", "ms", "tlen"], ascending=False)
        # for now just keep the first entry for each rbh
        df = df.drop_duplicates(subset=["rbh"])
        chromdf = df[["query", "target", "strand", "tstart", "tend"]]
        # sort by target, then tstart, then tend
        chromdf = chromdf.sort_values(by=["target", "tstart", "tend"])
        # save the chromdf as a tab delimited file
        chromdf.to_csv(output.chrom, sep="\t", index=False, header=False)

        # now we generate the protein fasta file
        with open(output.pep, "w") as o:
            for record in fasta.parse(input.pep):
                if record.id in chromdf["query"].values:
                    o.write(">{}\n{}\n".format(record.id, record.seq))

rule gzChrom:
    """
    just gzips the chrom file
    """
    input:
        chrom = config["tool"] + "/output/source_data/unannotated_genomes/{assemAnn}/{assemAnn}_annotated_with_{LG_name}.chrom",
        pep   = config["tool"] + "/output/source_data/unannotated_genomes/{assemAnn}/{assemAnn}_annotated_with_{LG_name}.pep"
    output:
        chrom = config["tool"] + "/output/source_data/unannotated_genomes/{assemAnn}/{assemAnn}_annotated_with_{LG_name}.chrom.gz",
        pep   = config["tool"] + "/output/source_data/unannotated_genomes/{assemAnn}/{assemAnn}_annotated_with_{LG_name}.pep.gz"
    threads: 1
    resources:
        mem_mb  = 1000, # shouldn't take a lot of RAM.
        time    = 5 # 5 minutes
    shell:
        """
        # first gzip the chrom file
        gzip < {input.chrom} > {output.chrom}
        # second gzip the pep file. This will be bigger.
        gzip < {input.pep} > {output.pep}
        """

rule generate_assembled_config_entry:
    """
    Print out a small piece of a yaml file specifically for ODP.
    These will be gathered and concatenated later.
    """
    input:
        unannot_genomes = config["unannotated_genome_chr_tsv"],
        genome          = config["tool"] + "/output/source_data/unannotated_genomes/{assemAnn}/{assemAnn}.chr.fasta.gz",
        protein = config["tool"] + "/output/source_data/unannotated_genomes/{assemAnn}/{assemAnn}_annotated_with_{LG_name}.pep.gz",
        chrom   = config["tool"] + "/output/source_data/unannotated_genomes/{assemAnn}/{assemAnn}_annotated_with_{LG_name}.chrom.gz"
    output:
        yaml   = config["tool"] + "/output/source_data/unannotated_genomes/{assemAnn}/{assemAnn}_annotated_with_{LG_name}.yaml.part",
    threads: 1
    resources:
        time    = 5,
        mem_mb  = 1000
    run:
        # load in the dataframe of the annotated genomes
        df = pd.read_csv(input.unannot_genomes, sep="\t")
        # strip leading and trailing whitespace from the column names because pandas can screw up sometimes
        df.columns = df.columns.str.strip()

        row = df.loc[df["Assembly Accession"] == wildcards.assemAnn]
        taxid = row["Organism Taxonomic ID"].values[0]

        # s is the output string.
        # h is headspace. Currently 2 spaces because the odp yaml files are compoased like so:
        #species:
        #  Sp_name:
        #    taxid:
        #    genus:
        #    species:
        #    assembly_accession:
        #    proteins:
        #    chrom:
        #    genome:
        #    minscaflen:
        #    ....

        # genus
        try:
            genus = row["Organism Name"].values[0].split(" ")[0]
        except:
            genus = row["Organism Name"].values[0]
        # species
        try:
            species = row["Organism Name"].values[0].split(" ")[1]
        except:
            species = "None"
        # cleanup the species name
        species = species.replace("sp.", "sp")
        assemAnn = wildcards.assemAnn

        # now we strip unwanted characters from the genus, species name, and assembly accession
        strip_these_chars = [" ", "_", "-"]
        for char in strip_these_chars:
            genus = genus.replace(char, "")
            species = species.replace(char, "")
            assemAnn = assemAnn.replace(char, "")

        minscaflen = 100000
        if "|" in [genus, species, taxid, wildcards.assemAnn]:
            raise ValueError("You cannot have a | character in the genus, species, taxid, or assembly accession.")
        # We currently cannot have "_" characters in the species name because of problems it causes downstream with snakemake.
        spstring = "{}{}-{}-{}".format(genus, species, taxid, assemAnn)
        # If there are more than 3 '-' characters in the string, then we have a problem.
        if spstring.count("-") > 3:
            raise ValueError("There are too many '-' characters in the string: {}".format(spstring))
        h = "  "
        s = ""
        s += h + "{}:\n".format(spstring)
        s += h + h + "assembly_accession: {}\n".format(wildcards.assemAnn)
        s += h + h + "taxid:              {}\n".format(taxid)
        s += h + h + "genus:              {}\n".format(genus)
        s += h + h + "species:            {}\n".format(species)
        s += h + h + "proteins:           {}\n".format(os.path.abspath(input.protein))
        s += h + h + "chrom:              {}\n".format(os.path.abspath(input.chrom))
        s += h + h + "genome:             {}\n".format(os.path.abspath(input.genome))
        s += h + h + "minscaflen:         {}\n".format(minscaflen)

        # write the string to the output
        with open(output.yaml, "w") as f:
            f.write(s)

rule collate_assembled_config_entries:
    """
    Concatenate all of the yaml files together into one big file.
    """
    input:
        yaml_parts = expand(config["tool"] + "/output/source_data/unannotated_genomes/{assemAnn}/{assemAnn}_annotated_with_{{LG_name}}.yaml.part",
                            assemAnn=config["assemAnn"])
    output:
        yaml = "NCBI_odp_db.unannotated.{LG_name}.yaml"
    threads: 1
    resources:
        mem_mb  = 1000,
        time    = 5
    shell:
        """
        echo "species:" > {output.yaml}
        cat {input.yaml_parts} >> {output.yaml}
        """

rule generate_species_list_for_timetree:
    """
    Generates a list of species that are in this database.
    One species per line.
    """
    input:
        yaml    = "NCBI_odp_db.unannotated.{LG_name}.yaml"
    output:
        sp_list = "NCBI_odp_sp_list.unannotated.{LG_name}.txt"
    threads: 1
    resources:
        time    = 5, # 5 minutes
        mem_mb  = 1000
    run:
        # open the yaml file into a dictionary
        with open(input.yaml, "r") as f:
            yaml_dict = yaml.load(f, Loader=yaml.FullLoader)
        # open the output file for writing
        with open(output.sp_list, "w") as f:
            # loop through the dictionary and write the species names to the file
            for sp in yaml_dict["species"]:
                f.write("{}\t{}\n".format(yaml_dict["species"][sp]["genus"], yaml_dict["species"][sp]["species"]))
