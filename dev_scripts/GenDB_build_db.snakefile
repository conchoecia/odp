"""
This takes the list of annotated and unannotated genomes and prepares a database from them for ODP.

"""

configfile: "config.yaml"

config["tool"] = "odp_build_db"

wildcard_constraints:
    taxid="[0-9]+",

import os

# first we must load in all of the files. Only do it once
rule all:
    input:
        config["tool"] + "/output/annotated_genomes.tsv",
        config["tool"] + "/output/unannotated_genomes.tsv"

# Purgatory
        ## safe mkdir for directories. DO NOT TOUCH THE DIRECTORY
        #create_directories_recursive_notouch(output.annotated_dir)
        #create_directories_recursive_notouch(output.unannotated_dir)

        ## make a directory for each assembly we want to download
        ## and make a yaml with its information. Do this by iterating through the dataframes
        ##
        ## Aug 4th 2023 - I realized that at some point there will need to be a control sequence
        ##  To recognize whether a genome has been updated from "unannotated" to "annotated".
        ##  This will especially be useful in cases where the unannotated genome is uploaded to NCBI.
        #for index, row in annot_df.iterrows():
        #    thisassembly = row["Assembly Accession"]
        #    outpath = os.path.join(output.annotated_dir, thisassembly)
        #    create_directories_recursive_notouch(outpath)
        ## now do the same as above for the unannotated genomes
        #for index, row in unann_df.iterrows():
        #    thisassembly = row["Assembly Accession"]
        #    outpath = os.path.join(output.unannotated_dir, thisassembly)
        #    create_directories_recursive_notouch(outpath)

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