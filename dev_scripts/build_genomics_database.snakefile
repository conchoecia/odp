"""
build a database of chromosome-scale genome assemblies using the NCBI website
"""
import numpy as np
import os
import pandas as pd

# import fasta parser from dependencies
snakefile_path = os.path.dirname(os.path.realpath(workflow.snakefile))
dependencies_path = os.path.join(snakefile_path, "../dependencies/fasta-parser")
sys.path.insert(1, dependencies_path)
import fasta

bin_path = os.path.join(snakefile_path, "../bin")

configfile: "config.yaml"

config["tool"] = "odp_database_builder"

wildcard_constraints:
    taxid="[0-9]+",

import os

def create_directories_recursive_notouch(path):
    """
    Unlike os.makedirs, this function will not touch a directory if it already exists.
    This is useful for snakemake because it will not re-run a rule if the output already exists.
    """
    parts = os.path.normpath(path).split(os.path.sep)
    current_path = ""

    for part in parts:
        current_path = os.path.join(current_path, part)
        if not os.path.exists(current_path):
            os.mkdir(current_path)
 
def get_assembly_datapacks(wildcards):
    checkpoint_results = list(checkpoints.split_into_annotated_and_unannotated.get(**wildcards).output)
    annotated_dir = checkpoint_results[0]  # this forces the checkpoint to be executed before we continue
    unannotated_dir = checkpoint_results[1]  # this forces the checkpoint to be executed before we continue

    # get the annotated assemblies
    REJECT = [".snakemake_timestamp"]

    assemAnn = [os.path.join(annotated_dir,  "{0}/{0}.yaml.part".format(assembly)) \
                for assembly in os.listdir(annotated_dir) \
                if assembly not in REJECT]
    return assemAnn
    
rule all:
    input:
        get_assembly_datapacks

#rule all:
#    input:
#        config["tool"] + "/output/source_data/"
#        #expand(config["tool"] + "/input/{taxid}.tsv",
#        #       taxid = config["taxids"]),
#        #expand(config["tool"] + "/input/selected_genomes_{taxid}.tsv",
#        #       taxid = config["taxids"])

rule install_datasets:
    output:
        datasets = os.path.join(bin_path, "datasets")
    threads: 1
    resources:
        mem_mb = 1000
    shell:
        """
        curl -o datasets 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets'
        mv datasets {output.datasets}
        chmod +x {output.datasets}
        """

rule install_dataformat:
    output:
        dataformat = os.path.join(bin_path, "dataformat")
    threads: 1
    resources:
        mem_mb = 1000
    shell:
        """
        curl -o dataformat 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/dataformat'
        mv dataformat {dataformat}
        chmod +x {dataformat}
        """

rule generate_assembled_config_entry:
    """
    Print out a small piece of a yaml file specifically for ODP.
    These will be gathered and concatenated later.
    """
    input:
        annotated_genomes = ancient(config["tool"] + "/output/annotated_genomes.tsv"),
        report            = ancient(config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.report.txt"),
        genome            = ancient(config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.fasta"),
        protein           = ancient(config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.pep"),
        chrom             = ancient(config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.chrom"),
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
        
        spstring = "{}{}{}".format(genus, species, taxid)
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