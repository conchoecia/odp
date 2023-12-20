"""
Date: 2023-12-19 This is an updated version of the original script that focuses on genomes identified as chromosome-scale by NCBI.

Build a database of chromosome-scale genome assemblies using the NCBI website

PREREQUISITES:
  - This script requires the ete3 toolkit to be installed. This can be done with conda:
    https://anaconda.org/conda-forge/ete3
  - You must also preload the taxonomy toolkit with the following commands:
    http://etetoolkit.org/docs/latest/tutorial/tutorial_ncbitaxonomy.html
    ```
    from ete3 import NCBITaxa
    ncbi = NCBITaxa()
    ncbi.update_taxonomy_database()
    ```
"""

from ete3 import NCBITaxa
import numpy as np
import os
import pandas as pd
from datetime import datetime

import yaml
import time

from Newick_to_common_ancestors import yaml_file_legal as yaml_file_legal

# import fasta parser from dependencies
snakefile_path = os.path.dirname(os.path.realpath(workflow.snakefile))
bin_path = os.path.join(snakefile_path, "../bin")

# in theory we can just remove this
#dependencies_path = os.path.join(snakefile_path, "../dependencies/fasta-parser")
#sys.path.insert(1, dependencies_path)
#import fasta

configfile: "config.yaml"

if "datetime" not in config:
    config["datetime"] = datetime.now().strftime('%Y%m%d%H%M')

config["tool"] = "odp_ncbi_genome_scraper"

wildcard_constraints:
    taxid="[0-9]+",

import os

def create_directories_recursive_notouch(path):
    """
    Unlike os.makedirs, this function will not touch a directory if it already exists.
    This is useful for snakemake because it will not re-run a rule if the output already exists.
    """
    parts = os.path.normpath(path).split(os.path.sep)
    # Determine whether to remove the last part of the path.
    # Basically we have to determine if the last part is intended to be a path or a file.
    # If it is a file, then we need to remove it.
    file_endings = [".txt", ".tsv", ".csv"]
    end_is_file = False
    for ending in file_endings:
        if parts[-1].endswith(ending):
            end_is_file = True
            break
    if end_is_file:
        parts = parts[:-1]

    current_path = ""
    for part in parts:
        current_path = os.path.join(current_path, part)
        if not os.path.exists(current_path):
            os.mkdir(current_path)
    # safe return code if done
    return 0

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
        expand(config["tool"] + "/output/annotated_genomes_chr_{datetime}.tsv",
                datetime = config["datetime"]),
        expand(config["tool"] + "/output/annotated_genomes_nonchr_{datetime}.tsv",
                datetime = config["datetime"]),
        expand(config["tool"] + "/output/unannotated_genomes_chr_{datetime}.tsv",
                datetime = config["datetime"]),
        expand(config["tool"] + "/output/unannotated_genomes_nonchr_{datetime}.tsv",
                datetime = config["datetime"]),

rule install_datasets:
    output:
        datasets = os.path.join(bin_path, "datasets")
    threads: 1
    resources:
        time   = 5, # 5 minutes
        mem_mb = 1000
    shell:
        """
        curl -o datasets 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/linux-amd64/datasets'
        mv datasets {output.datasets}
        chmod +x {output.datasets}
        """

rule install_dataformat:
    output:
        dataformat = os.path.join(bin_path, "dataformat")
    threads: 1
    resources:
        time   = 5, # 5 minutes
        mem_mb = 1000
    shell:
        """
        curl -o dataformat 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/dataformat'
        mv dataformat {output.dataformat}
        chmod +x {output.dataformat}
        """

rule download_json:
    input:
        datasets = os.path.join(bin_path, "datasets"),
    output:
        genome_report = config["tool"] + "/input/{taxid}.json"
    threads: 1
    resources:
        time   = 5, # 5 minutes
        mem_mb = 1000
    shell:
        """
        #{input.datasets} summary genome taxon --assembly-level chromosome --as-json-lines {wildcards.taxid} > {output.genome_report}
        {input.datasets} summary genome taxon --as-json-lines {wildcards.taxid} > {output.genome_report}
        """

all_fields     = ["accession",
                  "ani-best-ani-match-ani",
                  "ani-best-ani-match-assembly",
                  "ani-best-ani-match-assembly_coverage",
                  "ani-best-ani-match-category",
                  "ani-best-ani-match-organism",
                  "ani-best-ani-match-type_assembly_coverage",
                  "ani-best-match-status",
                  "ani-category",
                  "ani-check-status",
                  "ani-comment",
                  "ani-submitted-ani-match-ani",
                  "ani-submitted-ani-match-assembly",
                  "ani-submitted-ani-match-assembly_coverage",
                  "ani-submitted-ani-match-category",
                  "ani-submitted-ani-match-organism",
                  "ani-submitted-ani-match-type_assembly_coverage",
                  "ani-submitted-organism",
                  "ani-submitted-species",
                  "annotinfo-busco-complete",
                  "annotinfo-busco-duplicated",
                  "annotinfo-busco-fragmented",
                  "annotinfo-busco-lineage",
                  "annotinfo-busco-missing",
                  "annotinfo-busco-singlecopy",
                  "annotinfo-busco-totalcount",
                  "annotinfo-busco-ver",
                  "annotinfo-featcount-gene-non-coding",
                  "annotinfo-featcount-gene-other",
                  "annotinfo-featcount-gene-protein-coding",
                  "annotinfo-featcount-gene-pseudogene",
                  "annotinfo-featcount-gene-total",
                  "annotinfo-method",
                  "annotinfo-name",
                  "annotinfo-pipeline",
                  "annotinfo-provider",
                  "annotinfo-release-date",
                  "annotinfo-release-version",
                  "annotinfo-report-url",
                  "annotinfo-software-version",
                  "annotinfo-status",
                  "assminfo-assembly-method",
                  "assminfo-atypicalis-atypical",
                  "assminfo-atypicalwarnings",
                  "assminfo-bioproject",
                  "assminfo-bioproject-lineage-accession",
                  "assminfo-bioproject-lineage-parent-accession",
                  "assminfo-bioproject-lineage-parent-accessions",
                  "assminfo-bioproject-lineage-title",
                  "assminfo-biosample-accession",
                  "assminfo-biosample-attribute-name",
                  "assminfo-biosample-attribute-value",
                  "assminfo-biosample-bioproject-accession",
                  "assminfo-biosample-bioproject-parent-accession",
                  "assminfo-biosample-bioproject-parent-accessions",
                  "assminfo-biosample-bioproject-title",
                  "assminfo-biosample-description-comment",
                  "assminfo-biosample-description-organism-common-name",
                  "assminfo-biosample-description-organism-infraspecific-breed",
                  "assminfo-biosample-description-organism-infraspecific-cultivar",
                  "assminfo-biosample-description-organism-infraspecific-ecotype",
                  "assminfo-biosample-description-organism-infraspecific-isolate",
                  "assminfo-biosample-description-organism-infraspecific-sex",
                  "assminfo-biosample-description-organism-infraspecific-strain",
                  "assminfo-biosample-description-organism-name",
                  "assminfo-biosample-description-organism-pangolin",
                  "assminfo-biosample-description-organism-tax-id",
                  "assminfo-biosample-description-title",
                  "assminfo-biosample-ids-db",
                  "assminfo-biosample-ids-label",
                  "assminfo-biosample-ids-value",
                  "assminfo-biosample-last-updated",
                  "assminfo-biosample-models",
                  "assminfo-biosample-owner-contact-lab",
                  "assminfo-biosample-owner-name",
                  "assminfo-biosample-package",
                  "assminfo-biosample-publication-date",
                  "assminfo-biosample-status-status",
                  "assminfo-biosample-status-when",
                  "assminfo-biosample-submission-date",
                  "assminfo-blast-url",
                  "assminfo-description",
                  "assminfo-level",
                  "assminfo-linked-assm-accession",
                  "assminfo-linked-assm-type",
                  "assminfo-name",
                  "assminfo-notes",
                  "assminfo-paired-assm-accession",
                  "assminfo-paired-assm-changed",
                  "assminfo-paired-assm-manual-diff",
                  "assminfo-paired-assm-name",
                  "assminfo-paired-assm-only-genbank",
                  "assminfo-paired-assm-only-refseq",
                  "assminfo-paired-assm-status",
                  "assminfo-refseq-category",
                  "assminfo-release-date",
                  "assminfo-sequencing-tech",
                  "assminfo-status",
                  "assminfo-submitter",
                  "assminfo-suppression-reason",
                  "assminfo-synonym",
                  "assminfo-type",
                  "assmstats-contig-l50",
                  "assmstats-contig-n50",
                  "assmstats-gaps-between-scaffolds-count",
                  "assmstats-gc-count",
                  "assmstats-gc-percent",
                  "assmstats-genome-coverage",
                  "assmstats-number-of-component-sequences",
                  "assmstats-number-of-contigs",
                  "assmstats-number-of-organelles",
                  "assmstats-number-of-scaffolds",
                  "assmstats-scaffold-l50",
                  "assmstats-scaffold-n50",
                  "assmstats-total-number-of-chromosomes",
                  "assmstats-total-sequence-len",
                  "assmstats-total-ungapped-len",
                  "checkm-completeness",
                  "checkm-completeness-percentile",
                  "checkm-contamination",
                  "checkm-marker-set",
                  "checkm-marker-set-rank",
                  "checkm-species-tax-id",
                  "checkm-version",
                  "current-accession",
                  "organelle-assembly-name",
                  "organelle-bioproject-accessions",
                  "organelle-description",
                  "organelle-infraspecific-name",
                  "organelle-submitter",
                  "organelle-total-seq-length",
                  "organism-common-name",
                  "organism-infraspecific-breed",
                  "organism-infraspecific-cultivar",
                  "organism-infraspecific-ecotype",
                  "organism-infraspecific-isolate",
                  "organism-infraspecific-sex",
                  "organism-infraspecific-strain",
                  "organism-name",
                  "organism-pangolin",
                  "organism-tax-id",
                  "source_database",
                  "type_material-display_text",
                  "type_material-label",
                  "wgs-contigs-url",
                  "wgs-project-accession",
                  "wgs-url"]

fields_to_print = ["accession",
                  "annotinfo-busco-complete",
                  "annotinfo-busco-duplicated",
                  "annotinfo-busco-fragmented",
                  "annotinfo-busco-lineage",
                  "annotinfo-busco-missing",
                  "annotinfo-busco-singlecopy",
                  "annotinfo-busco-totalcount",
                  "annotinfo-busco-ver",
                  "annotinfo-featcount-gene-non-coding",
                  "annotinfo-featcount-gene-other",
                  "annotinfo-featcount-gene-protein-coding",
                  "annotinfo-featcount-gene-pseudogene",
                  "annotinfo-featcount-gene-total",
                  "annotinfo-method",
                  "annotinfo-name",
                  "annotinfo-pipeline",
                  "annotinfo-provider",
                  "annotinfo-release-date",
                  "annotinfo-release-version",
                  "annotinfo-report-url",
                  "annotinfo-software-version",
                  "annotinfo-status",
                  "assminfo-assembly-method",
                  "assminfo-atypicalis-atypical",
                  "assminfo-atypicalwarnings",
                  "assminfo-biosample-accession",
                  "assminfo-biosample-bioproject-accession",
                  "assminfo-biosample-bioproject-parent-accession",
                  "assminfo-biosample-bioproject-parent-accessions",
                  "assminfo-biosample-bioproject-title",
                  "assminfo-biosample-description-comment",
                  "assminfo-biosample-description-organism-common-name",
                  "assminfo-biosample-description-organism-infraspecific-breed",
                  "assminfo-biosample-description-organism-infraspecific-cultivar",
                  "assminfo-biosample-description-organism-infraspecific-ecotype",
                  "assminfo-biosample-description-organism-infraspecific-isolate",
                  "assminfo-biosample-description-organism-infraspecific-sex",
                  "assminfo-biosample-description-organism-infraspecific-strain",
                  "assminfo-biosample-description-organism-name",
                  "assminfo-biosample-description-organism-pangolin",
                  "assminfo-biosample-description-organism-tax-id",
                  "assminfo-biosample-description-title",
                  "assminfo-biosample-ids-db",
                  "assminfo-biosample-last-updated",
                  "assminfo-biosample-models",
                  "assminfo-biosample-owner-contact-lab",
                  "assminfo-biosample-owner-name",
                  "assminfo-biosample-package",
                  "assminfo-biosample-publication-date",
                  "assminfo-biosample-status-status",
                  "assminfo-biosample-status-when",
                  "assminfo-blast-url",
                  "assminfo-description",
                  "assminfo-level",
                  "assminfo-linked-assm-accession",
                  "assminfo-linked-assm-type",
                  "assminfo-name",
                  "assminfo-notes",
                  "assminfo-paired-assm-accession",
                  "assminfo-paired-assm-changed",
                  "assminfo-paired-assm-manual-diff",
                  "assminfo-paired-assm-name",
                  "assminfo-paired-assm-only-genbank",
                  "assminfo-paired-assm-only-refseq",
                  "assminfo-paired-assm-status",
                  "assminfo-refseq-category",
                  "assminfo-release-date",
                  "assminfo-sequencing-tech",
                  "assminfo-status",
                  "assminfo-submitter",
                  "assminfo-type",
                  "assmstats-contig-l50",
                  "assmstats-contig-n50",
                  "assmstats-gaps-between-scaffolds-count",
                  "assmstats-gc-count",
                  "assmstats-gc-percent",
                  "assmstats-genome-coverage",
                  "assmstats-number-of-component-sequences",
                  "assmstats-number-of-contigs",
                  "assmstats-number-of-organelles",
                  "assmstats-number-of-scaffolds",
                  "assmstats-scaffold-l50",
                  "assmstats-scaffold-n50",
                  "assmstats-total-number-of-chromosomes",
                  "assmstats-total-sequence-len",
                  "assmstats-total-ungapped-len",
                  "current-accession",
                  "organelle-assembly-name",
                  "organelle-bioproject-accessions",
                  "organelle-description",
                  "organelle-infraspecific-name",
                  "organelle-submitter",
                  "organelle-total-seq-length",
                  "organism-common-name",
                  "organism-infraspecific-breed",
                  "organism-infraspecific-isolate",
                  "organism-infraspecific-sex",
                  "organism-name",
                  "organism-pangolin",
                  "organism-tax-id",
                  "source_database",
                  "type_material-display_text",
                  "wgs-contigs-url",
                  "wgs-project-accession",
                  "wgs-url"]

rule format_json_to_tsv:
    input:
        genome_report = config["tool"] + "/input/{taxid}.json",
        dataformat = os.path.join(bin_path, "dataformat")
    output:
        report_tsv = temp(config["tool"] + "/input/{taxid}.tsv")
    params:
        fields = ",".join(fields_to_print),
        #fields = ",".join(all_fields)
    threads: 1
    resources:
        time   = 5, # 5 minutes
        mem_mb = 1000
    shell:
        """
        {input.dataformat} tsv genome --inputfile {input.genome_report} --fields {params.fields} > {output.report_tsv}
        """

def append_to_list(inputlist, to_append):
    """
    figure out if the thing to append is an int or an iterable
    """
    # int or numpy.int64
    if isinstance(to_append, int) or isinstance(to_append, np.int64):
        return inputlist + [to_append]
    else:
        return inputlist + list(to_append)

def prefer_representative_genomes(df):
    """
    Prefers the representative genomes if there are multiple rows with the same WGS URL

    Works similarly to prefer_SRA_over_others()
    """
    # sort by WGS URL, then by Assembly Refseq Category. If one is representative, keep that row
    gb = df.groupby("WGS URL")
    indices_to_keep = []
    for name, group in gb:
        if len(group) == 1:
            indices_to_keep = append_to_list(indices_to_keep, group.index[0])
        else:
            if "representative genome" in list(group["Assembly Refseq Category"].unique()):
                # keep the rows with this value
                indices_to_keep = append_to_list(indices_to_keep, group.index[group["Assembly Refseq Category"] == "representative genome"])
            else:
                # keep everything
                indices_to_keep = append_to_list(indices_to_keep, group.index )
    return df.loc[indices_to_keep]

def prefer_SRA_over_others(df):
    """
    Sometimes other filters for duplicate assemblies do not completely remove duplicate entries.
    This function will preferentially select the SRA accession over the others.
    The column to look in is "Assembly BioSample Sample Identifiers Database"
    Works similarly to prefer_representative_genomes()
    """
    gb = df.groupby("WGS URL")
    indices_to_keep = []
    for name, group in gb:
        # just keep this row
        if len(group) == 1:
            indices_to_keep = append_to_list(indices_to_keep, group.index[0])
        else:
            FILT_COLUMN = "Assembly BioSample Sample Identifiers Database"
            KEEP_THIS   = "SRA"
            if KEEP_THIS in list(group[FILT_COLUMN].unique()):
                # keep the rows with this value
                indices_to_keep = append_to_list(indices_to_keep, group.index[group[FILT_COLUMN] == KEEP_THIS])
            else:
                # keep everything
                indices_to_keep = append_to_list(indices_to_keep, group.index)
    return df.loc[indices_to_keep]

def prefer_assemblies_with_no_superseded(df):
    """
    prefer assemblies that have not been superseded by another assembly for the same species
    """
    gb = df.groupby("Organism Taxonomic ID")
    indices_to_keep = []
    for name, group in gb:
        # just keep this row
        if len(group) == 1:
            indices_to_keep = append_to_list(indices_to_keep, group.index[0])
        else:
            FILT_COLUMN = "Assembly Notes"
            PREFERENTIALLY_GET_RID_OF_THIS   = "superseded by newer assembly for species"

            if (PREFERENTIALLY_GET_RID_OF_THIS in list(group[FILT_COLUMN].unique())) and (len(group[FILT_COLUMN].unique()) > 1):
                # keep the rows with this value
                indices_to_keep = append_to_list(indices_to_keep, group.index[group[FILT_COLUMN] != PREFERENTIALLY_GET_RID_OF_THIS])
            else:
                # keep everything
                indices_to_keep = append_to_list(indices_to_keep, group.index)
    return df.loc[indices_to_keep]

def remove_specific_GCAs(df, filepath_of_GCAs):
    """
    We often decide post-hoc that we do not want to include certain assemblies.
    These assemblies have specific GCA accessions, and by removing them from the dataframe
       we can prevent them from being included in the final database.

    The file at filepath_of_GCAs should be a text file with one GCA accession per line.
       In this file, if the line starts with a comment character then we simply ignore that line.
       In the context of odp, that file will look like this:

       ```
       # These are assemblies that are malformed on NCBI, or are not chromosome-scale.
       # Do not use inline-comments for this file. Just use one assembly per line.
       GCA_021556685.1
       GCA_013368085.1
       GCA_017607455.1
       GCA_905250025.1
       ```
    """
    # first check that filepath_of_GCAs exists
    if not os.path.exists(filepath_of_GCAs):
        raise Exception("The file {} does not exist".format(filepath_of_GCAs))

    assemblies_to_ignore = set()
    with open(filepath_of_GCAs, "r") as f:
        for line in f:
            # remove leading and trailing whitespace
            line = line.strip()
            # ignore lines that start with a comment character
            if line.startswith("#"):
                continue
            # ignore empty lines
            if len(line) == 0:
                continue
            # add this line to the set of assemblies to ignore
            assemblies_to_ignore.add(line)

    # Remove rows that have values in assemblies_to_ignore values in the "Assembly Accession" column
    return df.loc[~df["Assembly Accession"].isin(assemblies_to_ignore)]

def get_best_contig_L50_assembly(df):
    """
    In great anticlimactic fashion we now pick the assembly with the lowest contig L50

    In almost all cases this assembly is the best chromosome-scale assembly
    """
    gb = df.groupby("Organism Taxonomic ID")
    indices_to_keep = []
    for name, group in gb:
        # just keep this row
        if len(group) == 1:
            indices_to_keep = append_to_list(indices_to_keep, group.index[0])
        else:
            # sort the group by  "Assembly Stats Contig L50", ascending
            group = group.sort_values(by="Assembly Stats Contig L50", ascending=True)
            # just get the top row since it has the lowest contig L50, and probably the best assembly
            indices_to_keep = append_to_list(indices_to_keep, group.index[0])
    return df.loc[indices_to_keep]

rule get_representative_genomes:
    """
    This rule is responsible for parsing the genome report and selecting the representative genomes
     for each species.

    Currently this does not support multiple genomes for one species.
    """
    input:
        report_tsv = config["tool"] + "/input/{taxid}.tsv",
        assembly_ignore_list = os.path.join(snakefile_path, "assembly_ignore_list.txt")
    output:
        representative_genomes = config["tool"] + "/input/selected_genomes_{taxid}.tsv"
    params:
        taxid_prefix = config["tool"] + "/input/taxid_info/"
    threads: 1
    resources:
        time  = 5, # 5 minutes
        mem_mb = 1000
    run:
        # load in the dataframe
        df = pd.read_csv(input.report_tsv, sep="\t")
        # strip leading and trailing whitespace from the column names because pandas can screw up sometimes
        df.columns = df.columns.str.strip()
        print("number of species is {}".format(len(df["Organism Taxonomic ID"].unique())))
        print("len of raw df is {}".format(len(df)))

        # remove assemblies that are not chromosome-scale
        df = remove_specific_GCAs(df, input.assembly_ignore_list)

        # Don't do this anymore, as the assemblies may be contigs that are whole chromosome-length
        # 2023-12-19 - Tested the script without this and the next step, filtering for Scaffold N50, also removes the same assemblies.
        #              So, just rely on the numbers instead of NCBI's analysis of contig level.
        ## remove the assemblies that are simply contigs using "Assembly Level"
        #df = df.loc[df["Assembly Level"] != "Contig"]
        #print("len of df after filtering contigs is {}".format(len(df)))

        # 2023-12-20: Not doing this anymore. I will already pick the best assembly per species, so this is redundant
        ## remove assemblies with a scaffold N50 of 50kb.
        ## This is rather generous, but I may raise the number later
        #df = df.loc[df["Assembly Stats Scaffold N50"] >= 50000]
        #print("len of df after filtering for N50 is {}".format(len(df)))

        # Filter out some duplicate assemblies, preferentially select the ones that are listed as representative genomes
        df = prefer_representative_genomes(df)
        print("len of df after keeping the representative sequences is {}".format(len(df)))

        # For each species, prefer SRA repository over others.
        df = prefer_SRA_over_others(df)
        print("len of df after preferring SRA is {}".format(len(df)))

        # For each taxonomic ID, prefer assemblies that have not been superseded.
        df = prefer_assemblies_with_no_superseded(df)
        print("len of df after preferring non-superseded is {}".format(len(df)))

        # For each taxonomic ID, get the assembly with the lowest contig L50. This is almost invariably the best assembly for this species.
        df = get_best_contig_L50_assembly(df)
        print("len of df after keeping the lowest contig L50 assembly is {}".format(len(df)))

        # 2023-12-20 - Not doing this anymore. I will do some extra filtering steps later if I need to.
        ## Get rid of assemblies with more than 5000 scaffolds unless they are chromosome-scale
        ## This is also generous - there is no reason a genome assembly should have even 5000 scaffolds now.
        #df = df.loc[(df["Assembly Stats Number of Scaffolds"] <= 5000) | (df["Assembly Level"] == "Chromosome")]
        #print("len of df after filtering out non-chr-scale assemblies with more than 5000 scaffolds {}".format(len(df)))

        # convert data type of df["Organism Taxonomic ID"] to int
        df["Organism Taxonomic ID"] = df["Organism Taxonomic ID"].astype(int)

        # make a new column called Lineage. Get the NCBI Taxa lineage from ete3 NCBITaxa
        ncbi = NCBITaxa()
        taxid_dict = {taxid: ";".join([str(x) for x in ncbi.get_lineage(taxid)]) for taxid in list(df["Organism Taxonomic ID"].unique())}
        df["Lineage"] = df["Organism Taxonomic ID"].map(taxid_dict)

        # save the dataframe to the output
        df.to_csv(output.representative_genomes, sep="\t", index=False)

# this checkpoint triggers re-evaluation of the DAG
checkpoint split_into_annotated_and_unannotated_and_chr_nonchr:
    """
    This rule is responsible for parsing the genome report and selecting the representative genomes
     for each species.

    Currently this does not support multiple genomes for one species.
    """
    input:
        representative_genomes = expand(config["tool"] + "/input/selected_genomes_{taxid}.tsv", taxid = config["taxids"])
    output:
        # DO NOT CHANGE THE ORDER OF THESE FILES. THE FUNCTION get_assemblies(wildcards) DEPENDS ON IT
        annotated_genomes_chr       =       config["tool"] + "/output/annotated_genomes_chr_{datetime}.tsv",
        annotated_genomes_nonchr    =       config["tool"] + "/output/annotated_genomes_nonchr_{datetime}.tsv",
        unannotated_genomes_chr     =       config["tool"] + "/output/unannotated_genomes_chr_{datetime}.tsv",
        unannotated_genomes_nonchr  =       config["tool"] + "/output/unannotated_genomes_nonchr_{datetime}.tsv"
    threads: 1
    resources:
        time   = 5, # 5 minutes
        mem_mb = 1000
    run:
        list_of_annotated_chr_dfs      = []
        list_of_annotated_nonchr_dfs   = []
        list_of_unannotated_chr_dfs    = []
        list_of_unannotated_nonchr_dfs = []
        # load in the dataframe
        for thisfile in input.representative_genomes:
            df = pd.read_csv(thisfile, sep="\t")
            # strip leading and trailing whitespace from the column names because pandas can screw up sometimes
            df.columns = df.columns.str.strip()

            # annotated genomes. Get "Annotation Release Date" is not Nan
            annot_df = df.loc[df["Annotation Release Date"].notna()]
            annot_chrom_df    = annot_df.loc[annot_df["Assembly Level"] == "Chromosome"]
            annot_notchrom_df = annot_df.loc[annot_df["Assembly Level"] != "Chromosome"]
            list_of_annotated_chr_dfs.append(   annot_chrom_df)
            list_of_annotated_nonchr_dfs.append(annot_notchrom_df)

            # unannotated genomes
            unann_df = df.loc[~df["Annotation Release Date"].notna()]
            unannot_chrom_df    = unann_df.loc[unann_df["Assembly Level"] == "Chromosome"]
            unannot_notchrom_df = unann_df.loc[unann_df["Assembly Level"] != "Chromosome"]
            list_of_unannotated_chr_dfs.append(   unannot_chrom_df)
            list_of_unannotated_nonchr_dfs.append(unannot_notchrom_df)

        # put all the annotated dataframes together
        annot_chr_df    = pd.concat(list_of_annotated_chr_dfs)
        annot_nonchr_df = pd.concat(list_of_annotated_nonchr_dfs)
        unann_chr_df    = pd.concat(list_of_unannotated_chr_dfs)
        unann_nonchr_df = pd.concat(list_of_unannotated_nonchr_dfs)

        # Remove duplicate rows that may have been picked up by nested taxids.
        #   For example, if we pick "Metazoa" and "Arthropoda", there will be many duplicate rows.
        annot_chr_df    = annot_chr_df.drop_duplicates()
        annot_nonchr_df = annot_nonchr_df.drop_duplicates()
        unann_chr_df    = unann_chr_df.drop_duplicates()
        unann_nonchr_df = unann_nonchr_df.drop_duplicates()

        # save the dataframe to the output
        annot_chr_df.to_csv(     output.annotated_genomes_chr,      sep="\t", index=False)
        annot_nonchr_df.to_csv(  output.annotated_genomes_nonchr,   sep="\t", index=False)
        unann_chr_df.to_csv(     output.unannotated_genomes_chr,    sep="\t", index=False)
        unann_nonchr_df.to_csv(  output.unannotated_genomes_nonchr, sep="\t", index=False)



