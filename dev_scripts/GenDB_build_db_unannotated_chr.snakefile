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

if "directory" in config:
    # ensure that the user also hasn't specified the accession tsvs
    if "accession_tsvs" in config:
        raise IOError("You cannot provide both a directory of tsv files ('directory') and specify the exact path of the TSVs ('accession_tsvs'). Read the config file.")
    # ensure that the directory exists
    if not os.path.isdir(config["directory"]):
        raise IOError("The directory of TSV files you provided does not exist. {}".format(config["directory"]))
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
    config["assemAnn"] = GenDB.determine_genome_accessions(config["unannotated_genome_chr_tsv"])
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
        raise IOError("You cannot provide both a directory of tsv files ('directory') and specify the exact path of the TSVs ('accession_tsvs'). Read the config file.")
    # we haven't implemented this yet. I haven't found a use case where I would want to specifially pick a file path rather than just get the most recent one.
    raise NotImplementedError("We haven't implemented this yet. I haven't found a use case where I would want to specifially pick a file path rather than just get the most recent one.")

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

wildcard_constraints:
    taxid="[0-9]+",

# first we must load in all of the files. Only do it once
rule all:
    input:
        expand(config["tool"] + "/output/source_data/unannotated_genomes/{assemAnn}/{assemAnn}.chr.fasta.gz", assemAnn=config["assemAnn"]),
        #LG_outfiles,
        expand(config["tool"] + "/output/source_data/unannotated_genomes/{assemAnn}/{assemAnn}_annotated_with_{LG_name}.chrom",
               assemAnn=config["assemAnn"], LG_name=LG_to_db_directory_dict.keys()),
        expand(config["tool"] + "/output/source_data/unannotated_genomes/{assemAnn}/{assemAnn}_annotated_with_{LG_name}.pep",
               assemAnn=config["assemAnn"], LG_name=LG_to_db_directory_dict.keys()),
        expand("NCBI_odp_db.unannotated.{LG_name}.yaml",
               LG_name=LG_to_db_directory_dict.keys()),
        expand("NCBI_odp_sp_list.unannotated.{LG_name}.txt",
               LG_name=LG_to_db_directory_dict.keys()),

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
        pep  = config["tool"] + "/input/LG_proteins/{LG_name}.fasta",
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
        chrom = config["tool"] + "/output/source_data/unannotated_genomes/{assemAnn}/{assemAnn}_annotated_with_{LG_name}.chrom",
        pep   = config["tool"] + "/output/source_data/unannotated_genomes/{assemAnn}/{assemAnn}_annotated_with_{LG_name}.pep"
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

rule download_unzip:
    """
    We have selected the unannotated genomes to download.
    For these genomes we need to find a way to annotate them.

    To specifically download the chromosome-scale scaffolds, there is a series of commands with the NCBI datasets tool.
      I found this set of instructions after opening a ticket on NCBI's github page: https://github.com/ncbi/datasets/issues/298
    """
    input:
        datasets = os.path.join(bin_path, "datasets")
    output:
        assembly = temp(config["tool"] + "/output/source_data/unannotated_genomes/{assemAnn}/ncbi_dataset.zip"),
        fasta    = temp(config["tool"] + "/output/source_data/unannotated_genomes/{assemAnn}/{assemAnn}.chr.fasta")
    retries: 3
    params:
        outdir   = config["tool"] + "/output/source_data/unannotated_genomes/{assemAnn}/",
        APIstring = "" if "API_key" not in locals() else "--api-key {}".format(locals()["API_key"])
    threads: 1
    resources:
        mem_mb = 2000, # 1 GB of RAM
        time   = 20  # 20 minutes.
    shell:
        """
        ## Wait a random amount of time up to 2 minutes to avoid overloading the NCBI servers.
        #SLEEPTIME=$((1 + RANDOM % 120))
        #echo "Sleeping for $SLEEPTIME seconds to avoid overloading the NCBI servers."
        #sleep $SLEEPTIME

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

        # make the final assembly fasta file from the individual chromosome's .fna files
        find {params.outdir} -name "*.fna" -exec cat {{}} \\; > {output.fasta}
        # Remove the files that we no longer need.
        find {params.outdir} -name "*.fna" -exec rm {{}} \\;
        """

#rule filter_fasta:
#    """
#    This filters the fasta to only contain the chromosome-scale scaffolds.
#
#    This rule will be very bespoke for the different sources of genome assemblies.
#
#    Sources of genome assemblies and features that define them:
#      - Sanger Darwin Tree of Life Project:
#        - jsonl: {"accession":"GCA_940337035.1","assemblyInfo":{"assemblyLevel":"Chromosome","assemblyName":"PGI_AGRIOTES_LIN_V1","assemblyStatus":"current","assemblyType":"haploid","bioprojectAccession":"PRJEB47908","bioprojectLineage":[{"bioprojects":[{"accession":"PRJEB47908","title":"Agriotes lineatus genome and annotation from the Pest Genomics Initiative."}]}],"biosample":{"accession":"SAMEA13407341","attributes":[{"name":"ENA first public","value":"2022-06-14"},{"name":"ENA last update","value":"2022-06-14"},{"name":"ENA-CHECKLIST","value":"ERC000011"},{"name":"External Id","value":"SAMEA13407341"},{"name":"INSDC center alias","value":"ROTHAMSTED RESEARCH"},{"name":"INSDC center name","value":"ROTHAMSTED RESEARCH"},{"name":"INSDC first public","value":"2022-06-14T00:18:37Z"},{"name":"INSDC last update","value":"2022-06-14T00:18:37Z"},{"name":"INSDC status","value":"public"},{"name":"Submitter Id","value":"Agriotes_lineatus_genome"},{"name":"collection_date","value":"2020"},{"name":"common name","value":"click beetle"},{"name":"geo_loc_name","value":"Canada"},{"name":"sample_name","value":"Agriotes_lineatus_genome"}],"description":{"comment":"Genome assembly of Agriotes lineatus. HiFi data assembled using Hifiasm. Dovetail genomics did the assembly and annotation processes (Maker). Single individual, unknown sex used for HiFi PacBio but with low coverage x5 and Omni-C Illumina reads. Sample taken: Canada","organism":{"organismName":"Agriotes lineatus","taxId":292458},"title":"Genome assembly of Agriotes lineatus from the Pest Genomics Initiative."},"lastUpdated":"2022-07-15T14:46:57.000","models":["Generic"],"owner":{"name":"EBI"},"package":"Generic.1.0","publicationDate":"2022-06-14T00:00:00.000","sampleIds":[{"db":"SRA","value":"ERS11009814"}],"status":{"status":"live","when":"2022-06-16T08:48:17.480"},"submissionDate":"2022-06-15T10:35:13.640"},"blastUrl":"https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastSearch&PROG_DEF=blastn&BLAST_SPEC=GDH_GCA_940337035.1","refseqCategory":"representative genome","releaseDate":"2022-06-24","submitter":"ROTHAMSTED RESEARCH"},"assemblyStats":{"contigL50":4905,"contigN50":201347,"gcCount":"1143871914","gcPercent":36.0,"genomeCoverage":"6.0x","numberOfComponentSequences":14154,"numberOfContigs":23837,"numberOfScaffolds":14154,"scaffoldL50":10,"scaffoldN50":95282374,"totalNumberOfChromosomes":10,"totalSequenceLength":"3159863714","totalUngappedLength":"3159019153"},"currentAccession":"GCA_940337035.1","organism":{"organismName":"Agriotes lineatus","taxId":292458},"sourceDatabase":"SOURCE_DATABASE_GENBANK","wgsInfo":{"masterWgsUrl":"https://www.ncbi.nlm.nih.gov/nuccore/CALNHX000000000.1","wgsContigsUrl":"https://www.ncbi.nlm.nih.gov/Traces/wgs/CALNHX01","wgsProjectAccession":"CALNHX01"}}
#        - The chromosome-scale scaffolds appear to start with the prefix "OW", and look like "OW679194.1"
#        - The whole string for the chromosome-scale scaffolds is like this: ">OW679194.1 Agriotes lineatus genome assembly, chromosome: 1"
#        - The unplaced scaffolds start with the prefix "CA"
#        - The whole string for the unplaced scaffolds look like this: ">CALNHX010000001.1 Agriotes lineatus genome assembly, contig: Scaffold_11__1_contigs__length_1480898, whole genome shotgun sequence"
#      - ?? source
#        - chromosome-scale strings look like this: ">CM057117.1 Ailuropoda melanoleuca isolate CPB_GP_2021 chromosome 10, whole genome shotgun sequence"
#        - scaffold-level strings look like this:   ">JAJSAN010000990.1 Ailuropoda melanoleuca isolate CPB_GP_2021 Scaffold_1001, whole genome shotgun sequence"
#        - jsonl: {"accession":"GCA_029963865.1","assemblyInfo":{"assemblyLevel":"Chromosome","assemblyMethod":"FALCON v. May-2020","assemblyName":"CPB_AME_v1","assemblyStatus":"current","assemblyType":"haploid","bioprojectAccession":"PRJNA784095","bioprojectLineage":[{"bioprojects":[{"accession":"PRJNA784095","title":"Ailuropoda melanoleuca isolate:CPB_GP_2021 Genome sequencing and assembly"}]}],"biosample":{"accession":"SAMN23470118","attributes":[{"name":"isolate","value":"CPB_GP_2021"},{"name":"age","value":"20"},{"name":"sex","value":"female"},{"name":"tissue","value":"blood"}],"bioprojects":[{"accession":"PRJNA784095"}],"description":{"organism":{"organismName":"Ailuropoda melanoleuca","taxId":9646},"title":"Model organism or animal sample from Ailuropoda melanoleuca"},"lastUpdated":"2023-05-13T00:40:21.613","models":["Model organism or animal"],"owner":{"contacts":[{}],"name":"Sichuan Key Laboratory of Conservation Biology on Endangered Wildlife"},"package":"Model.organism.animal.1.0","publicationDate":"2023-05-13T00:40:21.613","sampleIds":[{"label":"Sample name","value":"GPv1"},{"db":"SRA","value":"SRS11275548"}],"status":{"status":"live","when":"2023-05-13T00:40:21.613"},"submissionDate":"2021-11-27T11:56:02.827"},"blastUrl":"https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastSearch&PROG_DEF=blastn&BLAST_SPEC=GDH_GCA_029963865.1","releaseDate":"2023-05-12","sequencingTech":"PacBio","submitter":"Sichuan Key Laboratory of Conservation Biology on Endangered Wildlife"},"assemblyStats":{"contigL50":27,"contigN50":28556066,"gcCount":"1038989641","gcPercent":42.0,"genomeCoverage":"113.9x","numberOfComponentSequences":1335,"numberOfContigs":2111,"numberOfScaffolds":1335,"scaffoldL50":8,"scaffoldN50":134169173,"totalNumberOfChromosomes":21,"totalSequenceLength":"2478972487","totalUngappedLength":"2475793874"},"currentAccession":"GCA_029963865.1","organism":{"commonName":"giant panda","infraspecificNames":{"isolate":"CPB_GP_2021","sex":"female"},"organismName":"Ailuropoda melanoleuca","taxId":9646},"sourceDatabase":"SOURCE_DATABASE_GENBANK","wgsInfo":{"masterWgsUrl":"https://www.ncbi.nlm.nih.gov/nuccore/JAJSAN000000000.1","wgsContigsUrl":"https://www.ncbi.nlm.nih.gov/Traces/wgs/JAJSAN01","wgsProjectAccession":"JAJSAN01"}}
#    """
#    input:
#        fasta = config["tool"] + "/output/source_data/unannotated_genomes/{assemAnn}/{assemAnn}.fasta",
#    output:
#        fasta = config["tool"] + "/output/source_data/unannotated_genomes/{assemAnn}/{assemAnn}.filt.fasta",
#    threads: 1
#    resources:
#        mem_mb = 1000, # 1 GB of RAM
#        time   = 10  # 10 minutes.
#    run:
#        # open up the fasta file and print out the biggest 50 scaffolds
#        scaf_to_size = {}
#        for record in fasta.parse(input.fasta):
#            scaf_to_size[record.id] = len(record.seq)
#        # sort the scaffolds by size
#        sorted_scafs = sorted(scaf_to_size, key=scaf_to_size.get, reverse=True)
#        # print the top 50 scaffolds and their size
#        for i in range(50):
#            print(sorted_scafs[i], scaf_to_size[sorted_scafs[i]])
#        sys.exit()

rule gzip_fasta_file:
    """
    In this step zip the fasta file to conserve space.
    """
    input:
        genome = config["tool"] + "/output/source_data/unannotated_genomes/{assemAnn}/{assemAnn}.chr.fasta",
    output:
        genome = config["tool"] + "/output/source_data/unannotated_genomes/{assemAnn}/{assemAnn}.chr.fasta.gz"
    threads: 1
    resources:
        mem_mb = 1000, # 1 GB of RAM
        time   = 50   # Usually takes less than 10 minutes. Just do 50 for exceptional cases. Exceptional cases usually take 150 minutes.
    params:
        outdir   = config["tool"] + "/output/source_data/unannotated_genomes/{assemAnn}/",
    shell:
        """
        echo "Gzipping the fasta file."
        gzip {input.genome}
        # remove the .fna files again, in case the last step failed
        find {params.outdir} -name "*.fna" -exec rm {{}} \\;
        """

rule generate_assembled_config_entry:
    """
    Print out a small piece of a yaml file specifically for ODP.
    These will be gathered and concatenated later.
    """
    input:
        unannot_genomes = config["unannotated_genome_chr_tsv"],
        genome          = config["tool"] + "/output/source_data/unannotated_genomes/{assemAnn}/{assemAnn}.chr.fasta.gz",
        protein         = config["tool"] + "/output/source_data/unannotated_genomes/{assemAnn}/{assemAnn}_annotated_with_{LG_name}.pep",
        chrom           = config["tool"] + "/output/source_data/unannotated_genomes/{assemAnn}/{assemAnn}_annotated_with_{LG_name}.chrom",
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

        minscaflen = 100000
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
