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
        config["tool"] + "/output/annotated_genomes.tsv",
        config["tool"] + "/output/unannotated_genomes.tsv"

        # taxids
        expand(config["tool"] + "/input/selected_genomes_{taxid}.tsv", taxid = config["taxids"]),
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
        mem_mb = 1000
    shell:
        """
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
        fields = ",".join(fields_to_print)
    threads: 1 
    resources:
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
    prefers the representative genomes if there are multiple rows with the same WGS URL
    
    Works similarly to prefer_SRA_over_others()
    """
    # sort by WGS URL, then by Assembly Refseq Category. If one is representative, keep that row
    gb = df.groupby("WGS URL")
    indices_to_keep = []
    for name, group in gb:
        #if "Aplysia californica" in list(group["Organism Name"].unique()):
        #    print(group)
        #    print(group["Assembly Refseq Category"])
        # just keep this row
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
        report_tsv = config["tool"] + "/input/{taxid}.tsv"
    output:
        representative_genomes = config["tool"] + "/input/selected_genomes_{taxid}.tsv"
    threads: 1
    resources:
        mem_mb = 1000
    run:
        # load in the dataframe
        df = pd.read_csv(input.report_tsv, sep="\t") 
        # strip leading and trailing whitespace from the column names because pandas can screw up sometimes
        df.columns = df.columns.str.strip()
        print("number of species is {}".format(len(df["Organism Taxonomic ID"].unique())))
        print("len of raw df is {}".format(len(df)))

        # remove the assemblies that are simply contigs using "Assembly Level"
        df = df.loc[df["Assembly Level"] != "Contig"]
        print("len of df after filtering contigs is {}".format(len(df)))

        # remove assemblies with a scaffold N50 of 50kb.
        # This is rather generous, but I may raise the number later
        df = df.loc[df["Assembly Stats Scaffold N50"] >= 50000]
        print("len of df after filtering for N50 is {}".format(len(df)))

        # Filter out some duplicate assemblies, preferentially select the ones that are listed as representative genomes
        df = prefer_representative_genomes(df)
        print("len of df after keeping the representative sequences is {}".format(len(df)))
        
        # prefer SRA repository over others
        df = prefer_SRA_over_others(df)
        print("len of df after preferring SRA is {}".format(len(df)))

        # prefer assemblies that have not been superseded
        df = prefer_assemblies_with_no_superseded(df)
        print("len of df after preferring non-superseded is {}".format(len(df)))
        
        # get the assembly with the lowest contig L50. This is almost invariably the best assembly for this species
        df = get_best_contig_L50_assembly(df)
        print("len of df after keeping the lowest contig L50 assembly is {}".format(len(df)))
        
        # Get rid of assemblies with more than 5000 scaffolds unless they are chromosome-scale
        # This is also generous - there is no reason a genome assembly should have even 5000 scaffolds now.
        df = df.loc[(df["Assembly Stats Number of Scaffolds"] <= 5000) | (df["Assembly Level"] == "Chromosome")]
        print("len of df after filtering out non-chr-scale assemblies with more than 5000 scaffolds {}".format(len(df)))

        # save the dataframe to the output
        df.to_csv(output.representative_genomes, sep="\t", index=False)

# this checkpoint triggers re-evaluation of the DAG
checkpoint split_into_annotated_and_unannotated:
    """
    This rule is responsible for parsing the genome report and selecting the representative genomes
     for each species.
    
    Currently this does not support multiple genomes for one species.
    """
    input:
        representative_genomes = expand(config["tool"] + "/input/selected_genomes_{taxid}.tsv", taxid = config["taxids"])
    output:
        # DO NOT CHANGE THE ORDER OF THESE FILES. THE FUNCTION get_assemblies(wildcards) DEPENDS ON IT
        annotated_genomes   =       config["tool"] + "/output/annotated_genomes.tsv",
        unannotated_genomes =       config["tool"] + "/output/unannotated_genomes.tsv"
    threads: 1
    resources:
        mem_mb = 1000
    run:
        list_of_annotated_dfs   = []
        list_of_unannotated_dfs = []
        # load in the dataframe
        for thisfile in input.representative_genomes:
            df = pd.read_csv(thisfile, sep="\t") 
            # strip leading and trailing whitespace from the column names because pandas can screw up sometimes
            df.columns = df.columns.str.strip()

            # annotated genomes. Get "Annotation Release Date" is not Nan
            annot_df = df.loc[df["Annotation Release Date"].notna()]
            list_of_annotated_dfs.append(annot_df)
        
            # unannotated genomes
            unann_df = df.loc[~df["Annotation Release Date"].notna()]
            list_of_unannotated_dfs.append(unann_df)
        
        # put all the annotated dataframes together
        annot_df = pd.concat(list_of_annotated_dfs)
        unann_df = pd.concat(list_of_unannotated_dfs)

        # remove duplicate rows that may have been picked up by nested taxids
        annot_df = annot_df.drop_duplicates()
        unann_df = unann_df.drop_duplicates()

        # save the dataframe to the output
        annot_df.to_csv(output.annotated_genomes,   sep="\t", index=False)
        unann_df.to_csv(output.unannotated_genomes, sep="\t", index=False) 
        
        # safe mkdir for directories. DO NOT TOUCH THE DIRECTORY
        create_directories_recursive_notouch(output.annotated_dir)
        create_directories_recursive_notouch(output.unannotated_dir)

        # make a directory for each assembly we want to download
        # and make a yaml with its information. Do this by iterating through the dataframes
        #
        # Aug 4th 2023 - I realized that at some point there will need to be a control sequence
        #  To recognize whether a genome has been updated from "unannotated" to "annotated".
        #  This will especially be useful in cases where the unannotated genome is uploaded to NCBI.
        for index, row in annot_df.iterrows():
            thisassembly = row["Assembly Accession"]
            outpath = os.path.join(output.annotated_dir, thisassembly)
            create_directories_recursive_notouch(outpath)
        # now do the same as above for the unannotated genomes
        for index, row in unann_df.iterrows():
            thisassembly = row["Assembly Accession"]
            outpath = os.path.join(output.unannotated_dir, thisassembly)
            create_directories_recursive_notouch(outpath)

rule download_annotated_genomes:
    """
    We have selected the annotated genomes to download. These are the easiest to handle since we don't have to annotate them ourselves.
    """
    input:
        assembly_dir = ancient(config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}"),
        datasets = ancient(os.path.join(bin_path, "datasets"))
    output:
        assembly = temp(config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/ncbi_dataset.zip")
    params:
        outdir   = config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/"
    threads: 1
    resources:
        mem_mb = 1000
    shell:
        """
        echo "DOWNLOAD {wildcards.assemAnn}"
        #cd {params.outdir}
        #{input.datasets} download genome accession {wildcards.assemAnn} --include genome,protein,gff3,gtf
        """

rule unzip_annotated_genomes:
    """
    We just downloaded the genome data packet. Now unzip it, delete the zip file, and rename things
    """
    input:
        assembly = ancient(config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/ncbi_dataset.zip")
    output:
        genome  = config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.fasta",
        protein = config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.pep",
        gff     = config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.gff",
        gtf     = config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.gtf"
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
        find {params.outdir} -name "*.gtf" -exec mv {{}} {output.gtf} \;
        """

rule prep_chrom_file_from_NCBI:
    """
    This takes the output files from the NCBI database and puts them into the file format we need for odp, clink, et cetera
    """
    input:
        genome   = ancient(config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.fasta"),
        protein  = ancient(config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.pep"),
        gff      = ancient(config["tool"] + "/output/source_data/annotated_genomes/{assemAnn}/{assemAnn}.gff"),
        chromgen = ancient(os.path.join(snakefile_path, "..", "scripts", "NCBIgff2chrom.py"))
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