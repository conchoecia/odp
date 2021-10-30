#!/usr/bin/env python

from multiprocessing import Pool, cpu_count
import pandas as pd
import numpy as np
import copy
import sys

def permute_10k(df):
    observations = {x: 0 for x in range(1, 5000)}
    num_permutations = 10000
    scafs = list(df.columns)
    for i in range(num_permutations):
        for thisscaf in df.columns:
            df[thisscaf] = np.random.permutation(df[thisscaf].values)
        subbed = df.groupby(scafs).size().reset_index(name = "counts")["counts"].value_counts().to_dict()
        for key in subbed:
            observations[key] += 1
        if i % 10 == 0:
               print("  - Finished {}/{} ({:.2f}%) analyses.  ".format(
                   i, num_permutations,
                   (i/num_permutations)*100), end = "\r")
    return observations


def main():
    myfile = "/data/user/dschultz/synteny_4way_HcaSpongeCapsRho/synteny_analysis/MBH/COW_EMU_HCA_RESLi_mutual_best_hits.tsv"
    species_string = myfile.split("/")[-1].replace("_mutual_best_hits.tsv", "")

    df = pd.read_csv(myfile, sep = "\t", index_col = 0)

    all_species = [x for x in species_string.split("_")]
    scafs = ["{}_scaf".format(x) for x in all_species]
    df = df[scafs]

    num_analyses = int(1000000/10000)
    values = (df.copy() for i in range(num_analyses))
    #print(len(values))

    observations = {x: 0 for x in range(1, 5000)}
    with Pool(50) as pool:
        res = pool.map(permute_10k, values)
        for thisdict in res:
            for key in thisdict:
                observations[key] += thisdict[key]

    for i in range(1, 20):
        print("{}\t{}".format(i, observations[i]))

if __name__ == '__main__':
    main()
