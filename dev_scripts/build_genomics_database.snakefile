"""
build a database of chromosome-scale genome assemblies using the NCBI website
"""

import os

snakefile_path = os.path.dirname(os.path.realpath(workflow.snakefile))
bin_path = os.path.join(snakefile_path, "../bin")

configfile: "config.yaml"

config["tool"] = "odp_database_builder"

rule all:
    input:
        os.path.join(bin_path, "datasets"),
        os.path.join(bin_path, "dataformat")
        #expand(config["tool"] + "/assembly_list.txt")


rule install_datasets:
    output:
        datasets = os.path.join(bin_path, "datasets")
    threads: 1
    resources:
        mem_mb = 1000
        """
        curl -o datasets 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets'
        mv datasets {output}
        chmod +x {output}
        """

rule install_dataformat:
    output:
        datasets = os.path.join(bin_path, "dataformat")
    threads: 1
    resources:
        mem_mb = 1000
        """
        curl -o dataformat 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/dataformat'
        mv dataformat {output}
        chmod +x {output}
        """
