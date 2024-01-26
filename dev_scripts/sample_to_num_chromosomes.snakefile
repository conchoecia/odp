"""
This takes a config file and outputs a tsv file of
"{}\t{}".format(sampleid, number of chromosomes)
"""

configfile: "config.yaml"

# use the fasta modiule to get the number of chromosomes
# This block imports fasta-parser as fasta
import os
import sys
snakefile_path = os.path.dirname(os.path.realpath(workflow.snakefile))
dependencies_path = os.path.join(snakefile_path, "../dependencies/fasta-parser")
sys.path.insert(1, dependencies_path)
import fasta

rule all:
    input:
        "species_chrom_counts.tsv"

rule count_chromosomes:
    input:
        genome = lambda wildcards: config["species"][wildcards.sample]["genome"]
    output:
        "chrom_counts/{sample}.chrom_counts"
    threads: 1
    resources:
        mem_mb  = 200,
        runtime = 10
    shell:
        """
        # if the filename ends in gzip, use zcat | grep
        if [[ {input.genome} == *.gz ]]; then
            zcat {input.genome} | grep "^>" | wc -l > {output}
        else
            grep "^>" {input.genome} | wc -l > {output}
        fi
        """

rule cat_results:
    input:
        expand("chrom_counts/{sample}.chrom_counts",
                sample=config["species"].keys())
    output:
        outfile = "species_chrom_counts.tsv"
    threads: 1
    resources:
        mem_mb  = 200,
        runtime = 10
    run:
        outhandle = open(output.outfile, "w")
        print("sample\tchromosomes", file = outhandle)
        for sample in config["species"].keys():
            infile = "chrom_counts/{}.chrom_counts".format(sample)
            with open(infile, "r") as f:
                chroms = f.read().strip()
                print("{}\t{}".format(sample, chroms), file = outhandle)
        outhandle.close()
