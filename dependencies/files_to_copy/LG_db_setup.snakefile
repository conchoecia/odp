"""
This script builds HMMs out of a local directory of a .rbh file
 and of a directory called ./aligned/ with protein alignments

"""

import os
import pandas as pd

# find the rbh file in the current directory
files = os.listdir("./")
if not len([x for x in files if x.endswith("rbh")]) > 0:
    raise IOError("There is no rbh file in this directory. Did you delete it by accident? Open a github issue, please")
rbhfile = [x for x in files if x.endswith("rbh")][0]
# find the directory of aligned fasta files
if not os.path.exists("aligned"):
     raise IOError("There is no aligned directory here in this directory. Did you delete it by accident? Open a github issue, please")

# come up with the species in the rbh
testdf = pd.read_csv(rbhfile, sep = "\t", index_col = None)
rbh_species = [x.replace("_scaf", "") for x in testdf.columns
               if x.endswith("_scaf")]
species_string = "_".join(sorted(rbh_species))
rbh_entries = list(set(list(testdf["rbh"])))

config["tool"] = "{}".format(rbhfile.strip(".rbh"))

if len(rbh_species) < 2:
    raise IOError("There must be at least two species in the rbh file.")

rule all:
    input:
        "{}.hmm".format(config["tool"]),
        "{}.check".format(config["tool"])

rule make_hmm:
    input:
        aligned =  "aligned/{rbh}.fasta"
    output:
        hmm     = temp("hmms/{rbh}.hmm")
    threads: 1
    shell:
        """
        hmmbuild {output.hmm} {input.aligned}
        """

rule cat_hmm:
    input:
        hmms = expand("hmms/{rbh}.hmm", rbh = rbh_entries)
    output:
        hmm = "{}.hmm".format(config["tool"])
    threads: 1
    shell:
        """
        cat {input.hmms} > {output.hmm}
        rm -rf hmms/
        """

rule verify_hmm_complete:
    input:
        hmm = "{}.hmm".format(config["tool"])
    output:
        check = "{}.check".format(config["tool"])
    threads: 1
    run:
        unseen = [x for x in rbh_entries]
        with open(input.hmm, "r") as f:
            for line in f:
                line = line.strip()
                if line and line.startswith("NAME  "):
                    entry = line.replace("NAME  ", "").strip()
                    try:
                        unseen.remove(entry)
                    except:
                        pass
        if len(unseen) != 0:
            raise IOError("There were some entries in the rbh file that were not in the hmm file: {}".format(unseen))
        else:
            outhandle = open(output.check, "w")
            print("all entries present in the hmm", file = outhandle)
            outhandle.close()
