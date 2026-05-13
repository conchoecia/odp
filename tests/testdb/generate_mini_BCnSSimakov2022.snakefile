"""
This script uses the installed BCnSSimakov2022 database
 to generate a database that is 10% smaller than the original.
"""

import os
import pandas as pd
import sys

# The path to the original database is two directories up in LG_db/BCnSSimakov2022
# path to this very file, using the os module
snakefile_path = os.path.dirname(os.path.realpath(workflow.snakefile))
# import rbh_tools
source_path = os.path.join(snakefile_path, "../../source")
sys.path.insert(1, source_path)
import rbh_tools

original_rbh     = os.path.join(snakefile_path, "../../LG_db/BCnSSimakov2022/BCnSSimakov2022.rbh")
orig_aligned_dir = os.path.join(snakefile_path, "../../LG_db/BCnSSimakov2022/aligned")

rule all:
    input:
        "mini_LG_db/mBCnSSimakov2022/mBCnSSimakov2022.rbh",
        "mini_LG_db/mBCnSSimakov2022/mBCnSSimakov2022.hmm",
        "mini_LG_db/mBCnSSimakov2022/mBCnSSimakov2022.check"

rule determine_10percent_dataset:
    """
    Create a 10% subsample of the RBH dataset.
    """
    input:
        rbh = original_rbh
    output:
        sub = "mini_LG_db/mBCnSSimakov2022/mBCnSSimakov2022.rbh"
    params:
        sub_size = 0.1
    run:
        import pandas as pd
        df = rbh_tools.parse_rbh(input.rbh)
        df = df.sort_values("rbh").reset_index(drop=True)
        df_sub = df.groupby("gene_group")
        group_results = [group.head(int(len(group) * params.sub_size)) for _, group in df_sub]
        df_sub = pd.concat(group_results)
        df_sub.reset_index(drop=True, inplace=True)
        df_sub.to_csv(output.sub, sep="\t", index=False)

def get_rbh_entries():
    """
    Extract RBH entries dynamically after the checkpoint is complete.
    """
    checkpoint_output = checkpoints.copy_fasta_files.get().output[0]  # Use checkpoint, not rules

    # Get all FASTA files that were actually copied
    existing_files = os.listdir(checkpoint_output)

    # Extract RBH names from filenames (remove .fasta extension)
    rbh_entries = [f.replace(".fasta", "") for f in existing_files if f.endswith(".fasta")]

    return rbh_entries

checkpoint copy_fasta_files:
    """
    Copy the required fasta files based on the sampled RBH entries.
    """
    input:
        sub = "mini_LG_db/mBCnSSimakov2022/mBCnSSimakov2022.rbh"
    output:
        directory("mini_LG_db/mBCnSSimakov2022/aligned")  # Store copied files here
    params:
        outdir = "mini_LG_db/mBCnSSimakov2022/aligned"
    run:
        import shutil

        df = pd.read_csv(input.sub, sep="\t")
        rbh_entries = df["rbh"].unique().tolist()  # Adjust column name if needed

        os.makedirs(params.outdir, exist_ok=True)

        for rbh in rbh_entries:
            original_fasta = os.path.join(orig_aligned_dir, f"{rbh}.fasta")
            target_fasta   = os.path.join(params.outdir,    f"{rbh}.fasta")
            if os.path.exists(original_fasta):
                shutil.copy(original_fasta, target_fasta)
            else:
                print(f"Warning: {original_fasta} does not exist. Skipping.")

rule make_hmm:
    """
    Generate HMMs from aligned FASTA files.
    """
    input:
        aligned = "mini_LG_db/mBCnSSimakov2022/aligned/{rbh}.fasta"
    output:
        hmm = temp("mini_LG_db/mBCnSSimakov2022/hmms/{rbh}.hmm")
    params:
        hmmdir = "mini_LG_db/mBCnSSimakov2022/hmms"
    shell:
        """
        mkdir -p {params.hmmdir}
        hmmbuild {output.hmm} {input.aligned}
        """

rule cat_hmm:
    """
    Concatenate all HMM files into a single database.
    """
    input:
        hmms=lambda wildcards: expand("mini_LG_db/mBCnSSimakov2022/hmms/{rbh}.hmm",
                                      rbh=get_rbh_entries())
    output:
        hmm = "mini_LG_db/mBCnSSimakov2022/mBCnSSimakov2022.hmm"
    shell:
        """
        cat {input.hmms} > {output.hmm}
        rm -rf hmms/
        """

rule verify_hmm_complete:
    input:
        hmm = "mini_LG_db/mBCnSSimakov2022/mBCnSSimakov2022.hmm",
        sub = "mini_LG_db/mBCnSSimakov2022/mBCnSSimakov2022.rbh"
    output:
        check = "mini_LG_db/mBCnSSimakov2022/mBCnSSimakov2022.check"
    threads: 1
    run:
        # get what is in the rbh file
        df = rbh_tools.parse_rbh(input.sub)
        rbh_entries = df["rbh"].unique().tolist()
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

