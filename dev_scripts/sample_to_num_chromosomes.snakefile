"""
This takes a config file and outputs a tsv file of
"{}\t{}".format(sampleid, number of chromosomes)
"""

configfile: "config.yaml"

rule all:
    input:
        "species_chrom_counts.tsv"

rule count_chromosomes:
    input:
        genome = configfile["species"]["{sample}"]["genome"]
    output:
        "chrom_counts/{sample}.chrom_counts"
    shell:
        "grep '>' {input.genome} | wc -l > {output}"
