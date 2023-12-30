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

rule dlChrs:
    """
    We have selected the annotated genomes to download. We only want the chromosome-scale scaffolds.

    To specifically download the chromosome-scale scaffolds, there is a series of commands with the NCBI datasets tool.
      I found this set of instructions after opening a ticket on NCBI's github page: https://github.com/ncbi/datasets/issues/298
    """
    input:
        datasets = os.path.join(bin_path, "datasets")
    output:
        assembly = temp(config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/genomeDl/ncbi_dataset.zip"),
        readme   = temp(config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/genomeDl/README.md"),
        fetch    = temp(config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/genomeDl/ncbi_dataset/fetch.txt"),
        fasta    = temp(config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.chr.fasta"),
    retries: 3
    params:
        outdir   = config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/genomeDl/",
        rmdir    = config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/genomeDl/ncbi_dataset/data/{assemAnn}/",
        APIstring = "" if "API_key" not in locals() else "--api-key {}".format(locals()["API_key"])
    threads: 1
    resources:
        mem_mb = 500, # 500 MB of RAM - it never uses much
        time   = 20  # 20 minutes.
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
        {input.datasets} download genome accession {wildcards.assemAnn} --chromosomes all \
            {params.APIstring} --dehydrated || true

        # now we try to unzip it
        unzip -o ncbi_dataset.zip

        # rehydrate the dataset
        {input.datasets} rehydrate --directory . --match chr

        # go back to the original directory
        cd $RETURNHERE

        # First we remove any .fna files that have the assembly name in them. The files that have chr01.fna, chr02.fna, etc. are the ones we want.
        # See the error on github here: https://github.com/ncbi/datasets/issues/298
        find {params.outdir} -name "*genomic.fna" -exec rm {{}} \\;
        # We also should remove the unplaced scaffolds that are sometimes downloaded with the chromosome-scale scaffolds.
        find {params.outdir} -name "*unplaced.scaf.fna" -exec rm {{}} \\;

        # Now that we have removed the files that we don't want, we are (hopefully) left with just the chromosome-scale scaffolds.
        # make the final assembly fasta file from the individual chromosome's .fna files
        find {params.outdir} -name "*.fna" -exec cat {{}} \\; > {output.fasta}
        # Remove the files that we no longer need.
        find {params.outdir} -name "*.fna" -exec rm {{}} \\;

        # remove a directory we no longer need
        rm -rf {params.rmdir}
        """

rule gzip_fasta_file:
    """
    In this step zip the fasta file to conserve space.
    """
    input:
        genome = config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.chr.fasta",
    output:
        genome = config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.chr.fasta.gz"
    threads: 1
    resources:
        mem_mb = 1000, # 1 GB of RAM
        time   = 50   # Usually takes less than 10 minutes. Just do 50 for exceptional cases. Exceptional cases usually take 150 minutes.
    params:
        outdir   = config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/",
    shell:
        """
        echo "Gzipping the fasta file."
        gzip {input.genome}
        # remove the .fna files again, in case the last step failed
        find {params.outdir} -name "*.fna" -exec rm {{}} \\;
        """

def verChr_get_mem_mb(wildcards, attempt):
    """
    The amount of RAM needed for the script depends on the size of the input genome.
    """
    attemptdict = {1: 1000,
                   2: 2000,
                   3: 4000,
                   4: 8000,
                   5: 16000,
                   6: 32000,
                   7: 64000,
                   8: 128000}
    return attemptdict[attempt]

rule verChr:
    """
    This rule does some basic checks on the chromosome-scale genome assemblies to make sure that they are valid.
    Uses the fasta module included in the odp package to do this.
      - Checks that the fasta entries appear no more than once.
    """
    input:
        fasta = config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.chr.fasta.gz"
    output:
        report = config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.chr.fasta.gz.CHECKED"
    threads: 1
    retries: 8
    resources:
        mem_mb = verChr_get_mem_mb, # 1 GB of RAM
        time   = 5  # 5 minutes.
    run:
        seen_entries = set()
        for record in fasta.parse(input.fasta):
            if record.id in seen_entries:
                raise ValueError("The fasta entry {} appears more than once in the file {}.".format(record.id, input.fasta))
            else:
                seen_entries.add(record.id)
        # if we get to this point, then the fasta file is valid. We can write the report file.
        with open(output.report, "w") as f:
            f.write("The fasta file {} is valid.\n".format(input.fasta))
            f.write("It contains {} entries.\n".format(len(seen_entries)))
            f.write("The entries are: {}\n".format(seen_entries))
            f.write("None of the entires occur more than once.\n")

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
        protein  = temp(config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.pep"),
        gff      = temp(config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.gff")
    retries: 3
    params:
        outdir   = config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/pepDl/",
        APIstring = "" if "API_key" not in locals() else "--api-key {}".format(locals()["API_key"])
    threads: 1
    resources:
        mem_mb = 500, # Usually only uses 100MB of RAM
        time   = 5  # 5 minutes.
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
        report = config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.chr.fasta.gz.CHECKED", # we need this to verify that the fasta file is legal
        protein  = config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.pep",
        gff      = config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.gff",
        chromgen = os.path.join(snakefile_path, "..", "scripts", "NCBIgff2chrom.py")
    output:
        chrom  = temp(config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.chrFilt.chrom"),
        pep    = temp(config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.chrFilt.pep"),
        report = config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.chrFilt.report.txt"
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
        genome_checked    = config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.chr.fasta.gz.CHECKED", # we need this to verify that the fasta file is legal
        protein           = config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.chrFilt.pep.gz",
        chrom             = config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.chrFilt.chrom.gz",
        report            = config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.chrFilt.report.txt",
    output:
        yaml   = config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.yaml.part",
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

        if "|" in [genus, species, taxid, wildcards.assemAnn]:
            raise ValueError("You cannot have a | character in the genus, species, taxid, or assembly accession.")
        spstring = "{}{}|{}|{}".format(genus, species, taxid, wildcards.assemAnn)
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
        s += h + h + "minscafsize:         {}\n".format(str(int(minscaflen)))

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
        yaml = "NCBI_odp_db.annotated.chr.yaml"
    threads: 1
    resources:
        mem_mb  = 1000
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
        sp_list = "NCBI_odp_sp_list.annotated.chr.txt"
    threads: 1
    resources:
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