#!/usr/bin/env python

import itertools
import pandas as pd
import numpy as np
import sys
from collections import Counter

rbhfile = "/Users/darrin/Downloads/BFL_EMU_RES_reciprocal_best_hits.rbh"
rbhfile = "/Users/darrin/Downloads/COW_EMU_HCA_RESLi_reciprocal_best_hits.rbh"
rbhfile = "/Users/darrin/Downloads/EMU_RES_reciprocal_best_hits.rbh"


def permute_n_times(df, num_permutations):
    """
    Counts the number of observations of each combination of chromosomes.
    """
    df = df.copy()
    df = df[[x for x in df.columns if "_scaf" in str(x)]]

    species = [x.replace("_scaf","") for x in list(df.columns)]
    scafs = [np.array(df[x]) for x in df.columns]
    observations = {y: Counter() for y in itertools.product(*[np.unique(x) for x in scafs])}
    print("num combos: {}".format(len(observations)), file = sys.stderr)
    old_observations = Counter()
    #sys.exit()
    print(df, file = sys.stderr)
    print(scafs, file = sys.stderr)
    for i in range(num_permutations):
        for ii in range(len(scafs)):
            scafs[ii] = np.random.permutation(scafs[ii])
        counts = Counter()
        for ii in range(len(scafs[0])):
            counts[tuple([scafs[j][ii] for j in range(len(species))])] += 1
        for key in counts:
            observations[key][counts[key]] += 1
        if i % 10 == 0:
               print("  - Finished {}/{} ({:.2f}%) analyses.  ".format(
                   i, num_permutations,
                   (i/num_permutations)*100), end = "\r", file = sys.stderr)
        # now build the distribution fot the old observations
        for j in range(max(counts.values())+1):
            old_observations[j] += 1

    print("", file = sys.stderr)
    FDR_obs = {x: [] for x in observations}
    for key in observations:
        total = sum(observations[key].values())
        observations[key][0] = num_permutations - total
        runsum = 0
        # write a range in reverse
        for i in range(max(observations[key].keys()), -1, -1):
            runsum = observations[key][i] + runsum
            FDR_obs[key].append(runsum)
        FDR_obs[key] = [x/num_permutations for x in FDR_obs[key][::-1]]
        FDR_obs[key] = {i: FDR_obs[key][i] for i in range(len(FDR_obs[key]))}
            
        #printstring = "\t".join([str(observations[key][x]) for x in range(max(observations[key].keys())+1)])
        #print("{}\t{}".format(key, printstring))

    FDR_obs_old = []
    for i in range(max(old_observations.keys())+1):
        #print(i, old_observations[i]/old_observations[0])
        FDR_obs_old.append(old_observations[i]/old_observations[0])
    FDR_obs_old = {i: FDR_obs_old[i] for i in range(len(FDR_obs_old))}
    minval = 1/num_permutations
    for i in range(5001):
        if i not in FDR_obs_old:
            FDR_obs_old[i] = minval
    return species, FDR_obs, FDR_obs_old


df = pd.read_csv(rbhfile, sep="\t", header=0)
# groupby all the columns that have _scaf in them
gb = df.groupby([x for x in df.columns if "_scaf" in str(x)]).size().reset_index(name="num_genes_in_chr_group")
gb = gb.sort_values(by="num_genes_in_chr_group", ascending=False).reset_index(drop=True)

num_permutations = 100000
species, FDR_obs, FDR_obs_old = permute_n_times(df, num_permutations)
gb["lookup"] = gb.apply(lambda row : tuple([row["{}_scaf".format(x)] for x in species]), axis = 1)
gb["old_FDR"] = gb["num_genes_in_chr_group"].map(FDR_obs_old)
gb["new_FDR"] = gb.apply(lambda row : FDR_obs[row["lookup"]][row["num_genes_in_chr_group"]] if row["num_genes_in_chr_group"] in FDR_obs[row["lookup"]] else 1/num_permutations, axis = 1)
print(gb)
# gb to tsv
gb.to_csv("odp_nway_rbh_newshuffling.tsv", sep="\t", index=False)


#print(FDR_obs)
#print(FDR_obs_old)
#print(permute_output)

#        observations = permute_n_times(df, params.num_permutations)
#        with open(output.alpha, "w") as f:
#            for key in observations:
#                print("{}\t{}".format(key, observations[key]),
#                      file = f)
#
#        observations = {x: 0 for x in range(1, 5000)}
#        for thisfile in input.fdr_results:
#            with open(thisfile, "r") as f:
#                for line in f:
#                    line = line.strip()
#                    if line:
#                        print(line)
#                        fields = [int(x) for x in line.split()]
#                        observations[fields[0]] += fields[1]
#
#        resultsDF = pd.Series(observations).to_frame()
#        resultsDF.reset_index(inplace = True)
#        resultsDF.columns = ["Num_Genes_In_Chr_Group", "Num_Permutations"]
#        resultsDF["Total_Tests"] = params.num_permutations
#        resultsDF["alpha"] = resultsDF["Num_Permutations"]/params.num_permutations
#        print(resultsDF)
#        resultsDF.to_csv(output.alpha, sep="\t", index = False)