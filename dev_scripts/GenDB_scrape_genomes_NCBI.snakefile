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

# Some specific NCBI taxids cause problems with the NCBI datasets tool.
# This one, GCA_900186335.3, causes a parsing error: https://github.com/ncbi/datasets/issues/300
hardcoded_ignore_accessions = ["GCA_900186335.3",
                               "GCA_000002165.1", # This is the 2009 Celera Genomics mouse genome. It was found to be contaminated, and too large.
                               "GCF_000002265.2", # This is the Celera Genomics rat genome. It is currently suppressed.
                               ]
import datetime
from datetime import timedelta
from ete3 import NCBITaxa
import itertools
import numpy as np
import os
import pandas as pd
from datetime import datetime

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
        # report of the final dataset
        expand(config["tool"] + "/input/report_{taxid}.tsv",
                taxid = config["taxids"]),
        # report and plot of the filtered dataset
        expand(config["tool"] + "/input/report_history_filtered_{taxid}.tsv",
                taxid = config["taxids"]),
        expand(config["tool"] + "/input/report_history_filtered_{taxid}.pdf",
                taxid = config["taxids"]),
        # report and plot of the unfiltered dataset. This is more appropriate for looking at the growth of NCBI genomes.
        expand(config["tool"] + "/input/report_history_raw_{taxid}.tsv",
                taxid = config["taxids"]),
        expand(config["tool"] + "/input/report_history_raw_{taxid}.pdf",
                taxid = config["taxids"])

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

def get_best_row_for_each_assembly_accession(df):
    """
    Sometimes after running the datasets program, the same Assembly Accession will be listed multiple times.
    For the annotated and chromosome-scale genomes, we want every assembly possible, including multiple assemblies per species.

    So, we need to filter out the duplicate assembly rows.
    """
    gb = df.groupby("Assembly Accession")
    indices_to_keep = []
    double_indices = []
    for name, group in gb:
        if len(group) == 1:
            # just keep this row, this accession occurs only once
            indices_to_keep = append_to_list(indices_to_keep, group.index[0])
        else:
            double_indices = append_to_list(double_indices, group.index)
            # sort the group by  "Assembly Stats Contig L50", ascending
            group = group.sort_values(by="Assembly Stats Contig L50", ascending=True)
            # just get the top row since it has the lowest contig L50, and probably the best assembly
            indices_to_keep = append_to_list(indices_to_keep, group.index[0])
    # print out the tsv of the groups that have multiple rows per Assembly Accession. Save it as double_indices.tsv
    df.loc[double_indices].to_csv("double_indices.tsv", sep="\t", index=False)
    sys.exit()

def return_stats_string(df):
    """
    Takes in a dataframe and returns a string with some stats about the dataframe.
    This is useful for QCing the progression of the dataframe through the pipeline.

    The example string is like this:
    'There are {} annotated, chr-scale genomes, and {} assembly accessions, for {} unique taxids."
    """
    return "There are {} dataframe rows, and {} assembly accessions (genomes), for {} unique taxids.".format(
        len(df),
        len(df["Assembly Accession"].unique()),
        len(df["Organism Taxonomic ID"].unique()))

def prefer_SRA_over_others(df, groupby_col = "Assembly Accession"):
    """
    Sometimes other filters for duplicate assemblies do not completely remove duplicate entries.
    This function will preferentially select the SRA accession over the others.
    The column to look in is "Assembly BioSample Sample Identifiers Database"
    Works similarly to prefer_representative_genomes()
    """
    # 20231228 - changed this to Assembly Accession instead of WGS URL - sometimes there is no URL
    gb = df.groupby(groupby_col)
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

def prefer_refseq(df, groupby_col = "Organism Taxonomic ID",
                  groupby_col2 = "Assembly Stats Number of Contigs",
                  groupby_col3 = "Assembly Stats Total Sequence Length"):
    """
    For each taxonomic ID, prefer the NCBI RefSeq assembly over the others.
    To do this, we find potential duplicate assemblies by grouping by taxonomic ID, the number of contigs,
      and the total sequence length.

    The RefSeq entries will look identical to other assemblies, except that they will have a different annotation.
    """
    gb = df.groupby([groupby_col, groupby_col2, groupby_col3])
    indices_to_keep = []
    for name, group in gb:
        # just keep this row
        if len(group) == 1:
            indices_to_keep = append_to_list(indices_to_keep, group.index[0])
        else:
            #if there is a value "NCBI RefSeq" in the column "Annotation Provider", then keep only those rows
            FILT_COLUMN = "Annotation Provider"
            KEEP_THIS   = "NCBI RefSeq"
            if KEEP_THIS in list(group[FILT_COLUMN].unique()):
                # keep the rows with this value
                indices_to_keep = append_to_list(indices_to_keep, group.index[group[FILT_COLUMN] == KEEP_THIS])
            else:
                # keep everything
                indices_to_keep = append_to_list(indices_to_keep, group.index)
    return df.loc[list(set(indices_to_keep))]

def prefer_assemblies_with_no_superseded(df, groupby_col = "Assembly Accession"):
    """
    prefer assemblies that have not been superseded by another assembly for the same species
    """
    gb = df.groupby(groupby_col)
    indices_to_keep = []
    indices_for_groups = []
    for name, group in gb:
        # just keep this row
        if len(group) == 1:
            indices_to_keep = append_to_list(indices_to_keep, group.index[0])
        else:
            # keep the rows that have multiple assemblies per species
            indices_for_groups = append_to_list(indices_for_groups, group.index)
            FILT_COLUMN = "Assembly Notes"
            PREFERENTIALLY_GET_RID_OF_THIS   = "superseded by newer assembly for species"

            if (PREFERENTIALLY_GET_RID_OF_THIS in list(group[FILT_COLUMN].unique())) and (len(group[FILT_COLUMN].unique()) > 1):
                # keep the rows with this value
                indices_to_keep = append_to_list(indices_to_keep, group.index[group[FILT_COLUMN] != PREFERENTIALLY_GET_RID_OF_THIS])
            else:
                # keep everything
                indices_to_keep = append_to_list(indices_to_keep, group.index)
    return df.loc[list(set(indices_to_keep))]

def prefer_assembly_with_higher_N50(df, groupby_col = "Organism Taxonomic ID"):
    """
    Prefer assemblies that have a higher N50, be it a higher contig N50 or a higher scaffold N50.
    """
    gb = df.groupby(groupby_col)
    indices_to_keep = []
    indices_for_groups = []
    for name, group in gb:
        # just keep this row
        if len(group) == 1:
            indices_to_keep = append_to_list(indices_to_keep, group.index[0])
        else:
            # For each row get the highest N50 number, be it the contig N50 or the scaffold N50
            # Only keep the one with the higher N50.
            # If there is a tie, then keep both.
            group_highest_N50 = 0
            group_highest_N50_index = -1
            for index, row in group.iterrows():
                this_N50 = max(row["Assembly Stats Contig N50"], row["Assembly Stats Scaffold N50"])
                if this_N50 > group_highest_N50:
                    group_highest_N50 = this_N50
                    group_highest_N50_index = index
            indices_to_keep = append_to_list(indices_to_keep, group_highest_N50_index)
    return df.loc[list(set(indices_to_keep))]

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
    # First check if the input type is an iterable or a filepath.
    # If it is a filepath, then we should read in the file and get the set of GCAs to ignore.
    # If it is an iterable, then we should just use that iterable.

    assemblies_to_ignore = set()
    if isinstance(filepath_of_GCAs, str):
        # IN THIS CASE WE HAVE PASSED A STRING. THIS IS A FILEPATH, BUT WE NEED TO CHECK IT.
        # first check that filepath_of_GCAs exists
        if not os.path.exists(filepath_of_GCAs):
            raise Exception("The file {} does not exist".format(filepath_of_GCAs))

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
    else:
        # IN THIS CASE WE HAVE PASSED AN ITERABLE. USE THE ITERABLE.
        # verify that it is a set or list of strings
        acceptable_types = [set, list]
        if not any(isinstance(filepath_of_GCAs, x) for x in acceptable_types):
            raise Exception("The input type for filepath_of_GCAs is {}, but it should be a filepath or a set or list of strings.".format(type(filepath_of_GCAs)))
        assemblies_to_ignore = set(x for x in filepath_of_GCAs)

    # Remove rows that have values in assemblies_to_ignore values in the "Assembly Accession" column
    return df.loc[~df["Assembly Accession"].isin(assemblies_to_ignore)]

def get_best_contig_L50_assembly(df, groupby_col = "Assembly Accession"):
    """
    In great anticlimactic fashion we now pick the assembly with the lowest contig L50

    In almost all cases this assembly is the best chromosome-scale assembly
    """
    gb = df.groupby(groupby_col)
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

def legal_True_final_group_df(df) -> bool:
    """
    For the final dataframes for everything that is chromosome-scale, there should be the same number of rows
    as there are unique assembly accessions.

    Returns True if the dataframe is legal, raises an exception if it is not.
    """
    if len(df) != len(df["Assembly Accession"].unique()):
        raise Exception("There are {} rows in the dataframe, but {} unique assembly accessions.".format(
            len(df),
            len(df["Assembly Accession"].unique())))
    return True

def filter_raw_genome_df(df, hardcoded_ignore_accessions, suppress_text = False) -> pd.DataFrame:
    """
    This function filters out the raw list of genomes into a final list of genome assemblies that we want to use.
    The basic algorithm is this:
        1. Select the genomes that are chromosome-scale and are annotated.
          - For now we even select multiple individuals per species. This means that sometimes both haplotypes get used from one individual.
        2. Select the genomes that are chromosome-scale and are unannotated.
          - We also select multiple individuals per species here, and even both haplotypes if there are haplotype-resolved assemblies.
          - We also select species that also occur in the chromosome-scale, annotated list. This allows for more sampling.
        3. Select the genomes that are not chromosome-scale and are annotated.
          - We do not select genomes of species that already have chromosome-scale genomes.
          - We only select one individual per species here. The rules are detailed below, but it basically is just the assembly with the highest N50 for that species.
        4. Select the genomes that are not chromosome-scale and are not annotated.
          - Again, we do not select genomes of species that had chromosome-scale genomes in steps 1 and 2.
          - We only select one individual per species here.

    Returns a pandas dataframe of the filtered genomes.
    """
    fileout = sys.stderr
    if suppress_text:
        fileout = open(os.devnull, 'w')
    # remove rows that are absolute duplicates
    df = df.drop_duplicates()

    # Create new columns called "chrscale" and "annotated". Doing it this way avoids a SettingWithCopyWarning.
    df = df.assign(chrscale=False, annotated=False)

    print("Number of species is {}".format(len(df["Organism Taxonomic ID"].unique())), file = fileout)
    print("Len of raw df is {}".format(len(df)), file = fileout)

    print("", file = fileout)
    print("*** GETTING THE CHR-SCALE, ANNOTATED ASSEMBLIES ***", file = fileout)
    # First we get the assemblies that have annotations and are chromosome-scale
    df_annot_chr = df.loc[df["Annotation Release Date"].notna()]
    df_annot_chr = df_annot_chr.loc[df_annot_chr["Assembly Level"] == "Chromosome"]
    print("  - Getting the genomes that are annotated and chromosome-scale", file = fileout)
    print("    - {}".format(return_stats_string(df_annot_chr)), file = fileout)
    # Filter out some duplicate assemblies, preferentially select the ones that are listed as representative genomes
    print("  - Getting the rows that are SRA versions of each assembly accession.", file = fileout)
    df_annot_chr = prefer_SRA_over_others(df_annot_chr, groupby_col = "Assembly Accession")
    print("    - {}".format(return_stats_string(df_annot_chr)), file = fileout)
    # For each taxonomic ID, get the assembly with the lowest contig L50.
    print("  - Getting the rows that have the lowest contig L50 for each assembly accession.", file = fileout)
    print("    This doesn't really do anything because at this point, these rows will be duplicates.", file = fileout)
    df_annot_chr = get_best_contig_L50_assembly(df_annot_chr, groupby_col = "Assembly Accession")
    print("    - {}".format(return_stats_string(df_annot_chr)), file = fileout)
    df_annot_chr["chrscale"] = True
    df_annot_chr["annotated"] = True
    legal_True_final_group_df(df_annot_chr)

    print("", file = fileout)
    print("*** GETTING THE CHR-SCALE, unANNOTATED ASSEMBLIES ***", file = fileout)
    # first we get the assemblies that have no annotations and are chromosome-scale
    df_unannot_chr = df.loc[~df["Annotation Release Date"].notna()]
    df_unannot_chr = df_unannot_chr.loc[df_unannot_chr["Assembly Level"] == "Chromosome"]
    print("  - Getting the genomes that are unannotated and chromosome-scale", file = fileout)
    print("    - {}".format(return_stats_string(df_unannot_chr)), file = fileout)
    # Filter out some duplicate assemblies, preferentially select the ones that are listed as representative genomes
    print("  - Getting the rows that are the SRA versions of each assembly accession.", file = fileout)
    df_unannot_chr = prefer_SRA_over_others(df_unannot_chr, groupby_col = "Assembly Accession")
    print("    - {}".format(return_stats_string(df_unannot_chr)), file = fileout)
    # For each taxonomic ID, get the assembly with the lowest contig L50.
    print("  - Getting the rows that have the lowest contig L50 for each assembly accession.", file = fileout)
    print("    This doesn't really do anything because at this point, these rows will be duplicates.", file = fileout)
    df_unannot_chr = get_best_contig_L50_assembly(df_unannot_chr, groupby_col = "Assembly Accession")
    print("    - {}".format(return_stats_string(df_unannot_chr)), file = fileout)
    df_unannot_chr["chrscale"] = True
    df_unannot_chr["annotated"] = False
    legal_True_final_group_df(df_unannot_chr)

    # Generate a set of the species for which we have found chromosome-scale genomes already.
    # We don't need to find sub-chromosome-scale genomes for these species.
    # Everything is int already, so we don't need to cast.
    set_of_species_chr_scale = set(df_annot_chr["Organism Taxonomic ID"]).union(set(df_unannot_chr["Organism Taxonomic ID"]))

    print("", file = fileout)
    print("*** GETTING THE non-CHR-SCALE, ANNOTATED ASSEMBLIES ***", file = fileout)
    # first we get the assemblies that have annotations and are not chromosome-scale
    print("  - Getting the genomes that are annotated and not chromosome-scale", file = fileout)
    df_annot_nonchr = df.loc[df["Annotation Release Date"].notna()]
    df_annot_nonchr = df_annot_nonchr.loc[df_annot_nonchr["Assembly Level"] != "Chromosome"]
    print("    - {}".format(return_stats_string(df_annot_nonchr)), file = fileout)
    # Now we remove genomes of species that have already been found in the chromosome-scale datasets.
    print("  - Removing genomes for species that already have chromosome-scale genomes.", file = fileout)
    print("    The information we will learn from non-chromosome-scale genomes is redundant.", file = fileout)
    df_annot_nonchr = df_annot_nonchr.loc[~df_annot_nonchr["Organism Taxonomic ID"].isin(set_of_species_chr_scale)]
    print("    - {}".format(return_stats_string(df_annot_nonchr)), file = fileout)
    # Prefer the RefSeq versions of the genomes
    print("  - Filtering duplicate entries to prefer the assembly with the NCBI RefSeq annotation.", file = fileout)
    df_annot_nonchr = prefer_refseq(df_annot_nonchr)
    print("    - {}".format(return_stats_string(df_annot_nonchr)), file = fileout)
    # Prefer the SRA versions of the genomes
    print("  - Getting the rows that are the SRA versions of each assembly accession.", file = fileout)
    df_annot_nonchr = prefer_SRA_over_others(df_annot_nonchr, groupby_col = "Assembly Accession")
    print("    - {}".format(return_stats_string(df_annot_nonchr)), file = fileout)
    # filter out the duplicates now
    print("  - Getting the rows that have the lowest contig L50 for each assembly accession.", file = fileout)
    print("    This doesn't really do anything because at this point, these rows will be duplicates.", file = fileout)
    df_annot_nonchr = get_best_contig_L50_assembly(df_annot_nonchr, groupby_col = "Assembly Accession")
    print("    - {}".format(return_stats_string(df_annot_nonchr)), file = fileout)
    # prefer genomes with the higest N50, be it the scaffold or contig N50
    print("  - For each species, getting the assembly with the highest N50, be it contig or scaffold.", file = fileout)
    df_annot_nonchr = prefer_assembly_with_higher_N50(df_annot_nonchr)
    print("    - {}".format(return_stats_string(df_annot_nonchr)), file = fileout)
    df_annot_nonchr["chrscale"] = False
    df_annot_nonchr["annotated"] = True
    legal_True_final_group_df(df_annot_nonchr)


    print("", file = fileout)
    print("*** GETTING THE non-CHR-SCALE, nonANNOTATED ASSEMBLIES ***", file = fileout)
    # first we get the assemblies that do not annotations and are not chromosome-scale
    print("  - Getting the genomes that are unannotated and not chromosome-scale", file = fileout)
    df_unannot_nonchr = df.loc[~df["Annotation Release Date"].notna()]
    df_unannot_nonchr = df_unannot_nonchr.loc[df_unannot_nonchr["Assembly Level"] != "Chromosome"]
    print("    - {}".format(return_stats_string(df_unannot_nonchr)), file = fileout)
    # Now we remove genomes of species that have already been found in the chromosome-scale datasets.
    print("  - Removing genomes for species that already have chromosome-scale genomes.", file = fileout)
    print("    The information we will learn from non-chromosome-scale genomes is redundant.", file = fileout)
    df_unannot_nonchr = df_unannot_nonchr.loc[~df_unannot_nonchr["Organism Taxonomic ID"].isin(set_of_species_chr_scale)]
    print("    - {}".format(return_stats_string(df_unannot_nonchr)), file = fileout)
    # Prefer the SRA versions of the genomes
    print("  - Getting the rows that are the SRA versions of each assembly accession.", file = fileout)
    df_unannot_nonchr = prefer_SRA_over_others(df_unannot_nonchr, groupby_col = "Assembly Accession")
    print("    - {}".format(return_stats_string(df_unannot_nonchr)), file = fileout)
    # filter out the duplicates now
    print("  - Getting the rows that have the lowest contig L50 for each assembly accession.", file = fileout)
    print("    This doesn't really do anything because at this point, these rows will be duplicates.", file = fileout)
    df_unannot_nonchr = get_best_contig_L50_assembly(df_unannot_nonchr, groupby_col = "Assembly Accession")
    print("    - {}".format(return_stats_string(df_unannot_nonchr)), file = fileout)
    # prefer genomes with the higest N50, be it the scaffold or contig N50
    print("  - For each species, getting the assembly with the highest N50, be it contig or scaffold.", file = fileout)
    df_unannot_nonchr = prefer_assembly_with_higher_N50(df_unannot_nonchr)
    print("    - {}".format(return_stats_string(df_unannot_nonchr)), file = fileout)
    df_unannot_nonchr["chrscale"] = False
    df_unannot_nonchr["annotated"] = False
    legal_True_final_group_df(df_annot_nonchr)

    if suppress_text:
        fileout.close()

    # combine all of these into a new dataframe called df, line the original
    return pd.concat([df_annot_chr, df_unannot_chr, df_annot_nonchr, df_unannot_nonchr])

def load_and_cleanup_NCBI_datasets_tsv_df(tsv_filepath, ignore_list) -> pd.DataFrame:
    """
    This function is responsible for cleaning the dataframe that is output by the NCBI datasets program.
    A lot of times the fields are improperly formatted.
    """
    df = pd.read_csv(tsv_filepath, sep="\t", low_memory=False)
    # strip leading and trailing whitespace from the column names because pandas can screw up sometimes
    df.columns = df.columns.str.strip()

    # Remove assemblies that are not chromosome-scale.
    # These genomes are likely those that we manually checked and know that we don't wan't.
    df = remove_specific_GCAs(df, ignore_list)

    # remove duplicate rows
    df = df.drop_duplicates()

    # change these columns to integers. There should not be any missing values.
    cols_to_change_to_int = [
                             "Assembly Stats Contig L50",
                             "Assembly Stats Contig N50",
                             "Assembly Stats Number of Component Sequences",
                             "Assembly Stats Number of Contigs",
                             "Assembly Stats Total Sequence Length",
                             "Assembly Stats Total Ungapped Length",
                             "Organism Taxonomic ID",
                            ]
    for thiscol in cols_to_change_to_int:
        # check that the column exists, if not, tell the user
        if thiscol not in df.columns:
            raise Exception("The column {} is not in the dataframe.".format(thiscol))
        df[thiscol] = df[thiscol].astype(int)

    # change these columns to ints, and if there is NaN the value is 0
    cols_to_change_to_int = ["Annotation BUSCO Total Count",
                             "Annotation Count Gene Non-coding",
                             "Annotation Count Gene Other",
                             "Annotation Count Gene Protein-coding",
                             "Annotation Count Gene Pseudogene",
                             "Annotation Count Gene Total",
                             "Assembly Stats Number of Scaffolds",
                             "Assembly Stats Scaffold L50",
                             "Assembly Stats Scaffold N50"
                            ]
    for thiscol in cols_to_change_to_int:
        # check that the column exists, if not, tell the user
        if thiscol not in df.columns:
            raise Exception("The column {} is not in the dataframe.".format(thiscol))
        df[thiscol] = df[thiscol].fillna(0).astype(int)
    return df

def dataset_summary_table(df, assembly_release_date = "9999-99-99") -> pd.DataFrame:
    """
    This function takes in a dataframe and returns a summary table of the dataframe.

    The final dataframe looks like this, and has all the information needed for marginal tables:
            coltype chrscale annotated  num_species  num_genomes assembly_release_date_cutoff
            cell    False     False         4745         4745                   2023-12-29
            cell    False      True          851          851                   2023-12-29
            cell     True     False         2163         2837                   2023-12-29
            cell     True      True          648          789                   2023-12-29
        marginal    False       all         5248         5596                   2023-12-29
        marginal     True       all         2287         3626                   2023-12-29
        marginal      all     False         6908         7582                   2023-12-29
        marginal      all      True         1499         1640                   2023-12-29
           total      all       all         7535         9222                   2023-12-29
    """
    # if the assembly_release date is 9999-99-99, then return today's date
    if assembly_release_date == "9999-99-99":
        assembly_release_date = datetime.today().strftime('%Y-%m-%d')

    # print out a marginal table of the dataframes, the intersection of annotated, not annotated, chromosome-scale and not, plus the number of species in each category
    # the sort order of chrscale, then annotated
    #
    all_dfs = []
    summary_entries = []
    aggs = [["Organism Taxonomic ID", "nunique", "num_species"],
            ["Organism Taxonomic ID", "count", "num_genomes"]
           ]
    for thisagg in aggs:
        summary = df.groupby(["chrscale", "annotated"]).agg({thisagg[0]: thisagg[1]})
        # renmae the columns
        summary.columns = [thisagg[2]]
        summary_entries.append(summary)
    summary = pd.concat(summary_entries, axis=1).reset_index()
    # If any of these rows are missing from the dataframe, then add them with 0s
    # Specifically, the summary df must have values for all combinations of True/False for chrscale and annotated
    #  chrscale annotated  num_species  num_genomes
    #  False     False         4745         4745
    #  False      True          851          851
    #   True     False         2163         2837
    #   True      True          648          789
    colcombos = list(itertools.product([True, False], repeat=2))
    for thiscombo in colcombos:
        if (thiscombo[0], thiscombo[1]) not in list(zip(summary["chrscale"], summary["annotated"])):
            # make a dataframe of len 1 with the values of this combo
            tempdf = pd.DataFrame({"chrscale": thiscombo[0],
                                   "annotated": thiscombo[1],
                                   "num_species": 0,
                                   "num_genomes": 0},
                                  index=[0])
            # update the df
            summary = pd.concat([summary, tempdf]).reset_index(drop=True)
    summary = summary.assign(coltype="cell")
    all_dfs.append(summary)
    # DONE WITH THE CELL TABLE

    # NOW WE CALCULATE THE MARGINALS
    opposite = {"chrscale": "annotated", "annotated": "chrscale"}
    for marginal in ["chrscale", "annotated"]:
        marginal_entries = []
        for thisagg in aggs:
            summary = df.groupby([marginal]).agg({thisagg[0]: thisagg[1]})
            # rename the columns
            summary.columns = [thisagg[2]]
            marginal_entries.append(summary)
        summary = pd.concat(marginal_entries, axis=1).reset_index()
        # add the other column
        summary = summary.assign(**{opposite[marginal]: "all", "coltype": "marginal"})
        # At this point, every single time we should have a dataframe that looks like this:
        # When the marginal is chrscale, the summary df should look like this:
        #    chrscale  num_species  num_genomes annotated   coltype
        # 0     False         5269         5652       all  marginal
        # 1      True         2287         3633       all  marginal
        # To be explicit, the marginal column should have a "True" and "False row",
        #   the opposite column should say "all", every value in "coltype" should be "marginal"
        for thismarginal in [True, False]:
            if thismarginal not in list(summary[marginal]):
                # make a dataframe of len 1 with the values of this combo
                tempdf = pd.DataFrame({marginal: thismarginal,
                                       "num_species": 0,
                                       "num_genomes": 0,
                                       opposite[marginal]: "all",
                                       "coltype": "marginal"},
                                      index=[0])
                # update the df
                summary = pd.concat([summary, tempdf]).reset_index(drop=True)
        all_dfs.append(summary)

    # now just get the number of genomes and number of species for the whole dataframe
    # make a new dataframe with the number of genomes and number of species
    totdf = pd.DataFrame({"num_genomes": len(df["Assembly Accession"].unique()),
                          "num_species": len(df["Organism Taxonomic ID"].unique()),
                          "chrscale": "all",
                          "annotated": "all",
                          "coltype": "total"}, index=[0])
    all_dfs.append(totdf)
    # swap the columns so that coltype is first, then chrscale, then annotated, then num_species, then num_genomes
    all_dfs = [x[["coltype", "chrscale", "annotated", "num_species", "num_genomes"]] for x in all_dfs]

    # print out the concatenated dataframe
    finaldf = pd.concat(all_dfs).sort_values(by=["coltype", "chrscale", "annotated"]).reset_index(drop=True)
    # add a column of the assembly_release_date_cutoff
    finaldf = finaldf.assign(assembly_release_date_cutoff = assembly_release_date)
    return finaldf

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
        report                 = config["tool"] + "/input/report_{taxid}.tsv",
        representative_genomes = config["tool"] + "/input/selected_genomes_{taxid}.tsv"
    threads: 1
    resources:
        time  = 5, # 5 minutes
        mem_mb = 1000
    run:
        # first we get the list of things to ignore
        # remove things that are in assembly_ignore_list.txt
        ignore_list = []
        with open(input.assembly_ignore_list, "r") as f:
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
                ignore_list.append(line)
        # add hardcoded_ignore_accessions to the ignore_list
        ignore_list = set(ignore_list + hardcoded_ignore_accessions)

        # load in and clean up the dataframe
        df = load_and_cleanup_NCBI_datasets_tsv_df(input.report_tsv, ignore_list)

        # Perform filtering of the dataset to make sure that we have the best possible genomes.
        df = filter_raw_genome_df(df, hardcoded_ignore_accessions)

        # print out the dataset summary table
        summarydf = dataset_summary_table(df)
        summarydf.to_csv(output.report, sep="\t", index=False)

        # make a new column called Lineage. Get the NCBI Taxa lineage from ete3 NCBITaxa
        ncbi = NCBITaxa()
        taxid_dict = {taxid: ";".join([str(x) for x in ncbi.get_lineage(taxid)]) for taxid in list(df["Organism Taxonomic ID"].unique())}
        df["Lineage"] = df["Organism Taxonomic ID"].map(taxid_dict)

        # save the dataframe to the output
        df.to_csv(output.representative_genomes, sep="\t", index=False)

rule history_of_assemblies_filtered:
    """
    Takes in the genome report and outputs the stats of the dataframe for different points in time.
      - This rule is for the filtered version of the dataset, in which each species can occur in multiple categories.
      - This is better for looking at the datasets that will be used for whole-genome comparisons.
    Step backward 7 days in time until we run out of assemblies to consider.
    """
    input:
        report_tsv = config["tool"] + "/input/{taxid}.tsv",
    output:
        report     = config["tool"] + "/input/report_history_filtered_{taxid}.tsv",
    threads: 1
    resources:
        time  = 5, # 5 minutes
        mem_mb = 1000
    params:
        day_step = 7
    run:
        historical_view_dfs = []
        # load in and clean up the dataframe
        df = load_and_cleanup_NCBI_datasets_tsv_df(input.report_tsv, hardcoded_ignore_accessions)

        # go back in time 7 days at a time until we run out of assemblies to consider
        # go back until January 2000
        jan2000 = datetime.strptime("2000-01-01", '%Y-%m-%d')
        current_date = datetime.today().strftime('%Y-%m-%d')
        while current_date >= jan2000.strftime('%Y-%m-%d'):
            # print to sys.stderr in a progress-bar type configuration that just prints out the date
            print("   Filtering on or before date: {}".format(current_date), file=sys.stderr, end="\r")
            # filter the dataframe to only include assemblies that were released before seven_days_ago
            dftemp = df.loc[df["Assembly Release Date"] <= current_date]
            # annotations may have come at a later date
            dftemp.loc[dftemp["Annotation Release Date"] > current_date, "Annotation Release Date"] = np.nan
            # if there are no assemblies left, then we are done
            # print out the dataset summary table
            dftemp = filter_raw_genome_df(dftemp, [], suppress_text = True)
            summarydf = dataset_summary_table(dftemp, assembly_release_date = current_date)
            historical_view_dfs.append(summarydf)
            current_date = datetime.strptime(current_date, '%Y-%m-%d')
            current_date = (current_date - timedelta(days=params.day_step)).strftime('%Y-%m-%d')
        # print a newline to keep the progress bar on the screen/in the log
        print("   Filtering on or before date: {}".format(current_date), file=sys.stderr)
        # concatenate all the dataframes together
        concatdf = pd.concat(historical_view_dfs)
        # save to the output
        concatdf.to_csv(output.report, sep="\t", index=False)

rule assembly_report_plot_filtered:
    """
    Make a pdf of the filtered dataset assembly report.
    """
    input:
        report          = config["tool"] + "/input/report_history_filtered_{taxid}.tsv",
        plotting_script = os.path.join(snakefile_path, "../scripts/plot_NCBI_genomes_history.py")
    output:
        pdf             = config["tool"] + "/input/report_history_filtered_{taxid}.pdf",
    threads: 1
    resources:
        time  = 1, # 5 minutes
        mem_mb = 500
    shell:
        """
        python {input.plotting_script} -i {input.report} -o {output.pdf}
        """

rule history_of_assemblies_raw:
    """
    Takes in the genome report and outputs the stats of the dataframe for different points in time.
      - This rule is for the raw version of the dataset, in which each species can only occur in one category.
      - This is more accurate for looking at the absolute growth of the NCBI dataset over time.
    Step backward 7 days in time until we run out of assemblies to consider.
    """
    input:
        report_tsv = config["tool"] + "/input/{taxid}.tsv",
    output:
        report     = config["tool"] + "/input/report_history_raw_{taxid}.tsv",
    threads: 1
    resources:
        time  = 5, # 5 minutes
        mem_mb = 1000
    params:
        day_step = 7
    run:
        historical_view_dfs = []
        # load in and clean up the dataframe
        df = load_and_cleanup_NCBI_datasets_tsv_df(input.report_tsv, hardcoded_ignore_accessions)

        # go back in time 7 days at a time until we run out of assemblies to consider
        # go back until January 2000
        jan2000 = datetime.strptime("2000-01-01", '%Y-%m-%d')
        current_date = datetime.today().strftime('%Y-%m-%d')
        while current_date >= jan2000.strftime('%Y-%m-%d'):
            # print to sys.stderr in a progress-bar type configuration that just prints out the date
            print("   Filtering on or before date: {}".format(current_date), file=sys.stderr, end="\r")
            # filter the dataframe to only include assemblies that were released before seven_days_ago
            dftemp = df.loc[df["Assembly Release Date"] <= current_date]
            # if there are no assemblies left, then we are done
            # print out the dataset summary table
            # sort by Assembly Accession, then Annotation release date, preferring the ones with annotations first
            dftemp = dftemp.sort_values(by=["Assembly Accession", "Annotation Release Date"], ascending=[True, False])
            # drop duplicates, keeping the first one
            dftemp = dftemp.drop_duplicates(subset=["Assembly Accession"], keep="first")
            # add the columns called "chrscale" and "annotated", set them to False
            dftemp = dftemp.assign(chrscale=False, annotated=False)
            # for "Annotation Release Date", make the cells after current date NaN
            dftemp.loc[dftemp["Annotation Release Date"] > current_date, "Annotation Release Date"] = np.nan
            # if Annotation Release Date is NaN, then set "annotated" to False, otherwise True. We already handled False, just do True
            dftemp.loc[~dftemp["Annotation Release Date"].isna(), "annotated"] = True
            # if Assembly level == "Chromosome", then set "chrscale" to True
            dftemp.loc[dftemp["Assembly Level"] == "Chromosome", "chrscale"] = True
            summarydf = dataset_summary_table(dftemp, assembly_release_date = current_date)
            historical_view_dfs.append(summarydf)
            current_date = datetime.strptime(current_date, '%Y-%m-%d')
            current_date = (current_date - timedelta(days=params.day_step)).strftime('%Y-%m-%d')
        # print a newline to keep the progress bar on the screen/in the log
        print("   Filtering on or before date: {}".format(current_date), file=sys.stderr)
        # concatenate all the dataframes together
        concatdf = pd.concat(historical_view_dfs)
        # save to the output
        concatdf.to_csv(output.report, sep="\t", index=False)

rule assembly_report_plot_raw:
    """
    Make a pdf of the raw dataset assembly report.
    """
    input:
        report          = config["tool"] + "/input/report_history_raw_{taxid}.tsv",
        plotting_script = os.path.join(snakefile_path, "../scripts/plot_NCBI_genomes_history.py")
    output:
        pdf             = config["tool"] + "/input/report_history_raw_{taxid}.pdf",
    threads: 1
    resources:
        time  = 1, # 5 minutes
        mem_mb = 500
    shell:
        """
        python {input.plotting_script} -i {input.report} -o {output.pdf}
        """


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

            # annotated genomes.
            # use the "annotated" and "chrscale" columns to determine what goes into what
            annot_chrom_df    = df.loc[df["annotated"] == True].loc[df["chrscale"] == True]
            annot_notchrom_df = df.loc[df["annotated"] == True].loc[df["chrscale"] == False]
            list_of_annotated_chr_dfs.append(      annot_chrom_df)
            list_of_annotated_nonchr_dfs.append(annot_notchrom_df)

            # unannotated genomes.
            # use the "annotated" and "chrscale" columns to determine what goes into what
            unannot_chrom_df    = df.loc[df["annotated"] == False].loc[df["chrscale"] == True]
            unannot_notchrom_df = df.loc[df["annotated"] == False].loc[df["chrscale"] == False]
            list_of_unannotated_chr_dfs.append(      unannot_chrom_df)
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