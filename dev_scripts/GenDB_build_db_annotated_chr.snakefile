"""
This takes the list of annotated and unannotated genomes and prepares a database from them for ODP.

The program requires that either a directory of the annotated and unannotated genome lists be provided.

Otherwise the user has to specify specific paths to those tsv files.
"""

import os
import pandas as pd
import sys
from datetime import datetime
import GenDB
import yaml

#if "datetime" not in config:
#    config["datetime"] = datetime.now().strftime('%Y%m%d%H%M')

# figure out where bin is because we need to use some outside tools
snakefile_path = os.path.dirname(os.path.realpath(workflow.snakefile))
bin_path = os.path.join(snakefile_path, "../bin")

#if "API_key" not in locals():
#    API_key = ""

configfile: "config.yaml"
config["tool"] = "odp_ncbi_genome_db"
# Do some logic to see if the user has procided enough information for us to analyse the genomes
if ("directory" not in config) and ("accession_tsvs" not in config):
    raise IOerror("You must provide either a directory of the annotated and unannotated genome lists, or a list of the paths to those tsv files. Read the config file.")

if "directory" in config:
    # ensure that the user also hasn't specified the accession tsvs
    if "accession_tsvs" in config:
        raise IOerror("You cannot provide both a directory of tsv files ('directory') and specify the exact path of the TSVs ('accession_tsvs'). Read the config file.")
    # ensure that the directory exists
    if not os.path.isdir(config["directory"]):
        raise IOerror("The directory of TSV files you provided does not exist. {}".format(config["directory"]))
    # get the paths to the tsv files
    latest_accesions = GenDB.return_latest_accession_tsvs(config["directory"])
    # This is the output of the above function:
    #return_dict = {"annotated_chr":      os.path.join(directory_path, most_recent_annotated_chr_file),
    #               "annotated_nonchr":   os.path.join(directory_path, most_recent_annotated_nonchr_file),
    #               "unannotated_chr":    os.path.join(directory_path, most_recent_unannotated_chr_file),
    #               "unannotated_nonchr": os.path.join(directory_path, most_recent_unannotated_nonchr_file)
    #               }
    config["annotated_genome_chr_tsv"]       = latest_accesions["annotated_chr"]
    config["annotated_genome_nonchr_tsv"]    = latest_accesions["annotated_nonchr"]
    config["unannotated_genome_chr_tsv"]     = latest_accesions["unannotated_chr"]
    config["unannotated_genome_nonchr_tsv"]  = latest_accesions["unannotated_nonchr"]

    # now add the entries to the config file so we can download them or not
    config["assemAnn"] = GenDB.determine_genome_accessions(config["annotated_genome_chr_tsv"])

    # get the list of GCAs to ignore in case we need to remove any
    ignore_list_path = os.path.join(snakefile_path, "assembly_ignore_list.txt")
    with open(ignore_list_path, "r") as f:
        for line in f:
            line = line.strip()
            if line:
                if line in config["assemAnn"]:
                    config["assemAnn"].remove(line)

elif "accession_tsvs" in config:
    # ensure that the user also hasn't specified the directory
    if "directory" in config:
        raise IOerror("You cannot provide both a directory of tsv files ('directory') and specify the exact path of the TSVs ('accession_tsvs'). Read the config file.")
    # we haven't implemented this yet. I haven't found a use case where I would want to specifially pick a file path rather than just get the most recent one.
    raise NotImplementedError("We haven't implemented this yet. I haven't found a use case where I would want to specifially pick a file path rather than just get the most recent one.")

#onstart:
#    # NOTE: onstart doesn't execute if the workflow is run with the -n flag
#    # If we need to build a big database, we likely will need to use an API key
#    if config["require_API"] in [True, "true", "True", "TRUE", 1]:
#        # now that we're sure that we want to use an API key, make sure that the user has not saved
#        #  one to their file. This is not a secure way to do this. Instead we will prompt the user to type it in.
#        if "API_key" in config:
#            raise ValueError("You have specified that you want to use an API key, but you have saved one to your config file. This is not secure. Please remove the API key from your config file and run the program again. You will be prompted to enter your API key.")
#        else:
#            # prompt the user to enter their API key
#            global API_key
#            API_key = input("Please enter your NCBI API key then press enter: ")
#    else:
#        config["require_API"] = False

wildcard_constraints:
    taxid="[0-9]+",

# first we must load in all of the files. Only do it once
rule all:
    input:
        #expand(config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.yaml.part", assemAnn=config["assemAnn"]),
        "NCBI_odp_db.yaml",
        "NCBI_odp_sp_list.txt"

rule download_annotated_genomes:
    """
    We have selected the annotated genomes to download. These are the easiest to handle since we don't have to annotate them ourselves.
    """
    input:
        datasets = os.path.join(bin_path, "datasets")
    output:
        assembly = temp(config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/ncbi_dataset.zip")
    params:
        outdir   = config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/",
        APIstring = "" if "API_key" not in locals() else "--api-key {}".format(locals()["API_key"])
    threads: 1
    resources:
        mem_mb = 1000
    shell:
        """
        cd {params.outdir}
        {input.datasets} download genome accession {wildcards.assemAnn} {params.APIstring} --include genome,protein,gff3,gtf
        echo "Sleping for 5 minutes to avoid overloading the NCBI servers."
        sleep 300
        """

rule unzip_annotated_genomes:
    """
    We just downloaded the genome data packet. Now unzip it, delete the zip file, and rename things
    """
    input:
        assembly = config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/ncbi_dataset.zip"
    output:
        genome   = config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.fasta",
        protein  = config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.pep",
        gff      = config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.gff",
        #gtf      = config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.gtf"
    threads: 1
    resources:
        mem_mb = 1000
    params:
        outdir   = config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/"
    shell:
        """
        TMPDIR=`pwd`
        cd {params.outdir}
        unzip -o ncbi_dataset.zip
        cd $TMPDIR
        find {params.outdir} -name "*.fna" -exec mv {{}} {output.genome} \;
        find {params.outdir} -name "*.faa" -exec mv {{}} {output.protein} \;
        find {params.outdir} -name "*.gff" -exec mv {{}} {output.gff} \;
        """

rule prep_chrom_file_from_NCBI:
    """
    This takes the output files from the NCBI database and puts them into the file format we need for odp, clink, et cetera
    """
    input:
        genome   = config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.fasta",
        protein  = config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.pep",
        gff      = config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.gff",
        chromgen = os.path.join(snakefile_path, "..", "scripts", "NCBIgff2chrom.py")
    output:
        chrom   = config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.chrom",
        report  = config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.report.txt"
    threads: 1
    resources:
        mem_mb  = 1000
    params:
        prefix  = config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}"
    shell:
        """
        python {input.chromgen} -f {input.genome} -p {input.protein} -g {input.gff} -o {params.prefix}
        """

rule generate_assembled_config_entry:
    """
    Print out a small piece of a yaml file specifically for ODP.
    These will be gathered and concatenated later.
    """
    input:
        annotated_genomes = config["annotated_genome_tsv"],
        report            = config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.report.txt",
        genome            = config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.fasta",
        protein           = config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.pep",
        chrom             = config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.chrom",
    output:
        yaml   = config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.yaml.part",
    threads: 1
    resources:
        mem_mb  = 1000
    run:
        # load in the dataframe of the annotated genomes
        df = pd.read_csv(input.annotated_genomes, sep="\t")
        # strip leading and trailing whitespace from the column names because pandas can screw up sometimes
        df.columns = df.columns.str.strip()

        # read in the report into pandas if we can. Strip all the excess whitespace since the columns are formatted with whitespace
        reportdf = pd.read_csv(input.report, comment = "#", delim_whitespace=True)
        # get only the scaffolds that comprise more than 0.5% of the annotated proteins
        filtdf = reportdf.loc[reportdf["percent_of_proteins"] >= 0.05]
        minscaflen = filtdf["scaflen"].min() - 1000

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

        spstring = "{}{}{}".format(genus, species, taxid)
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
        yaml = "NCBI_odp_db.yaml"
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
        yaml    = "NCBI_odp_db.yaml"
    output:
        sp_list = "NCBI_odp_sp_list.txt"
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
