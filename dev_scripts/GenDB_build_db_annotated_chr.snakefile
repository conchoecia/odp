"""
This takes the list of annotated and unannotated genomes and prepares a database from them for ODP.

The program requires that either a directory of the annotated and unannotated genome lists be provided.

Otherwise the user has to specify specific paths to those tsv files.
"""

import pandas as pd
from datetime import datetime
import GenDB
import yaml

# This block imports fasta-parser as fasta
import os
import sys
snakefile_path = os.path.dirname(os.path.realpath(workflow.snakefile))
dependencies_path = os.path.join(snakefile_path, "../dependencies/fasta-parser")
sys.path.insert(1, dependencies_path)
import fasta

# figure out where bin is because we need to use some outside tools
bin_path = os.path.join(snakefile_path, "../bin")

configfile: "config.yaml"
config["tool"] = "odp_ncbi_genome_db"
config = GenDB.opening_logic_GenDB_build_db(config, chr_scale = True, annotated = True)

wildcard_constraints:
    taxid="[0-9]+",

# first we must load in all of the files. Only do it once
rule all:
    input:
        expand(config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.chr.fasta.gz",
            assemAnn = config["assemAnn"]),
        expand(config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.chrFilt.pep.gz",
            assemAnn = config["assemAnn"]),
        expand(config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.chrFilt.chrom.gz",
            assemAnn = config["assemAnn"]),
        "NCBI_odp_db.annotated.chr.yaml",
        "NCBI_odp_sp_list.annotated.chr.txt"

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
    We have selected the annotated genomes to download. We only want the chromosome-scale scaffolds.

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
       fasta   = temp(ensure(config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.chr.fasta", non_empty=True)),
       allscaf = ensure(config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.scaffold_df.all.tsv", non_empty=True),
       chrscaf = ensure(config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.scaffold_df.chr.tsv", non_empty=True)
    retries: 3
    params:
        outdir   = config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/",
    threads: 1
    group: "dlgz"
    resources:
        mem_mb = dlChrs_get_mem_mb, # the amount of RAM needed depends on the size of the input genome. Just scale UP.
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
        genome = config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.chr.fasta",
    output:
        genome = config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.chr.fasta.gz"
    threads: 1
    group: "dlgz"
    resources:
        mem_mb = 1000, # 1 GB of RAM
        time   = lambda wildcards: GenDB.gzip_get_time(config["assemAnn_to_scaflen"][wildcards.assemAnn])
    params:
        outdir   = config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/",
    shell:
        """
        echo "Gzipping the fasta file."
        gzip < {input.genome} > {output.genome}
        """

rule dlPepGff:
    """
    Because the other rule downloads only the chromosome-scale scaffolds, we need to download
      the pep file and gff file for this entry.

    Same structure as the previous download task.
    """
    input:
        datasets = os.path.join(bin_path, "datasets")
    output:
        readme   = temp(config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/pepDl/README.md"),
        assembly = temp(config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/pepDl/{assemAnn}.pepAndGff.zip"),
        protein  = temp(ensure(config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.pep", non_empty=True)),
        gff      = temp(ensure(config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.gff", non_empty=True))
    retries: 3
    params:
        outdir   = config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/pepDl/",
        APIstring = "" if "API_key" not in locals() else "--api-key {}".format(locals()["API_key"])
    threads: 1
    resources:
        mem_mb = 500, # Usually only uses 100MB of RAM
        time   = 5,  # 5 minutes.
        download_slots = 1
    shell:
        """
        # Wait a random amount of time up to 10 seconds to space out requests to the NCBI servers.
        SLEEPTIME=$((1 + RANDOM % 10))
        echo "Sleeping for $SLEEPTIME seconds to avoid overloading the NCBI servers."
        sleep $SLEEPTIME

        # Save the current directory
        # Function to download the file
        RETURNHERE=$(pwd)
        cd {params.outdir}
        {input.datasets} download genome accession {wildcards.assemAnn} \
            {params.APIstring} \
            --include protein,gff3,gtf \
            --filename {wildcards.assemAnn}.pepAndGff.zip || true

        # now we try to unzip it
        unzip -o {wildcards.assemAnn}.pepAndGff.zip

        # go back to the original directory
        cd $RETURNHERE

        # move the gff and faa files to the correct location
        find {params.outdir} -name "*.faa" -exec mv {{}} {output.protein} \\;
        find {params.outdir} -name "*.gff" -exec mv {{}} {output.gff} \\;
        # remove the gtf file if it exists
        find {params.outdir} -name "*.gtf" -exec rm {{}} \\;
        """

def prep_chrom_get_mem_mb(wildcards, attempt):
    """
    The amount of RAM needed for the script depends on the size of the input genome, the number of proteins, and the gff size.
    """
    attemptdict = {1: 4000,
                   2: 8000,
                   3: 16000,
                   4: 32000,
                   5: 64000,
                   6: 128000,
                   7: 256000,
                   8: 512000}
    return attemptdict[attempt]

def prep_chrom_get_time(wildcards, attempt):
    """
    The amount of minutes needed varies depending on the input size.
    """
    attemptdict = {1: 16,
                   2: 32,
                   3: 64,
                   4: 128,
                   5: 256,
                   6: 512,
                   7: 1024,
                   8: 2048}
    return attemptdict[attempt]

rule prep_chrom_file_from_NCBI:
    """
    This takes the output files from the NCBI database and puts them into the file format we need for odp, clink, et cetera

    This works with gzipped fasta files.
    """
    input:
        genome   = config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.chr.fasta.gz",
        protein  = config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.pep",
        gff      = config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.gff",
        chromgen = os.path.join(snakefile_path, "..", "scripts", "NCBIgff2chrom.py")
    output:
        chrom  = temp(ensure(config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.chrFilt.chrom", non_empty=True)),
        pep    = temp(ensure(config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.chrFilt.pep",   non_empty=True)),
        report = ensure(config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.chrFilt.report.txt", non_empty=True)
    threads: 1
    retries: 7
    resources:
        mem_mb  = prep_chrom_get_mem_mb, # shouldn't take much RAM, 231228 - I have seen mostly 200 MB or less. Sometimes it blows up to multiple GB.
        time    = prep_chrom_get_time    # Most of these end by 5 minutes, but occassionally they take longer.
    params:
        prefix  = config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.chrFilt"
    shell:
        """
        rm -f {output.chrom} {output.pep} {output.report} 2>/dev/null
        python {input.chromgen} \
            -f {input.genome} \
            -p {input.protein} \
            -g {input.gff} \
            --union \
            -o {params.prefix}
        """

rule gzPepGff:
    """
    This gzips the pep and gff files to conserve disk space.
    """
    input:
        chrom  = config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.chrFilt.chrom",
        pep    = config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.chrFilt.pep",
    output:
        chrom  = config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.chrFilt.chrom.gz",
        pep    = config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.chrFilt.pep.gz",
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
        annotated_genomes = config["annotated_genome_chr_tsv"],
        genome            = config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.chr.fasta.gz",
        protein           = config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.chrFilt.pep.gz",
        chrom             = config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.chrFilt.chrom.gz",
        report            = config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.chrFilt.report.txt",
    output:
        yaml   = ensure(config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.yaml.part", non_empty=True),
    threads: 1
    resources:
        mem_mb  = 1000,
        time    = 5
    run:
        # load in the dataframe of the annotated genomes
        df = pd.read_csv(input.annotated_genomes, sep="\t")
        # strip leading and trailing whitespace from the column names because pandas can screw up sometimes
        df.columns = df.columns.str.strip()

        # read in the report into pandas if we can. Strip all the excess whitespace since the columns are formatted with whitespace
        reportdf = pd.read_csv(input.report, comment = "#", delim_whitespace=True)
        minscaflen = reportdf["scaflen"].min() - 1000

        row = df.loc[df["Assembly Accession"] == wildcards.assemAnn]
        taxid = int(row["Organism Taxonomic ID"].values[0])

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
        strip_these_chars = [" ", "_", "-",]
        for char in strip_these_chars:
            genus = genus.replace(char, "")
            species = species.replace(char, "")
            assemAnn = assemAnn.replace(char, "")


        if "|" in [genus, species, taxid, wildcards.assemAnn]:
            raise ValueError("You cannot have a | character in the genus, species, taxid, or assembly accession.")
        # We currently cannot have "_" characters in the species name because of problems it causes downstream with snakemake.
        spstring = "{}{}-{}-{}".format(genus, species, taxid, assemAnn)
        # if there are more than three '-' characters in the string, then we have a problem.
        if spstring.count("-") > 3:
            raise ValueError("There are too many '-' characters in the string: {}".format(spstring))
        h = "  "
        s = ""
        s += h + "{}:\n".format(spstring)
        s += h + h + "assembly_accession: {}\n".format(wildcards.assemAnn)
        s += h + h + "taxid:              {}\n".format(str(int(taxid)))
        s += h + h + "genus:              {}\n".format(genus)
        s += h + h + "species:            {}\n".format(species)
        s += h + h + "proteins:           {}\n".format(os.path.abspath(input.protein))
        s += h + h + "chrom:              {}\n".format(os.path.abspath(input.chrom))
        s += h + h + "genome:             {}\n".format(os.path.abspath(input.genome))
        s += h + h + "minscafsize:        {}\n".format(str(int(minscaflen)))

        # write the string to the output
        with open(output.yaml, "w") as f:
            f.write(s)

rule collate_assembled_config_entries:
    """
    Concatenate all of the yaml files together into one big file.
    """
    input:
        yaml_parts = expand(config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.yaml.part",
                            assemAnn=config["assemAnn"])
    output:
        yaml = ensure("NCBI_odp_db.annotated.chr.yaml", non_empty=True)
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
        yaml    = "NCBI_odp_db.annotated.chr.yaml"
    output:
        sp_list = ensure("NCBI_odp_sp_list.annotated.chr.txt", non_empty = True)
    threads: 1
    resources:
        mem_mb  = 1000,
        time    = 5
    run:
        # open the yaml file into a dictionary
        with open(input.yaml, "r") as f:
            yaml_dict = yaml.load(f, Loader=yaml.FullLoader)
        # open the output file for writing
        with open(output.sp_list, "w") as f:
            # loop through the dictionary and write the species names to the file
            for sp in yaml_dict["species"]:
                f.write("{}\t{}\n".format(yaml_dict["species"][sp]["genus"], yaml_dict["species"][sp]["species"]))