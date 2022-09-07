#!/usr/bin/env python

"""
This code handles the coloring of the orthlogs in a rbh file.

There are a few ways to do this. The first way is to just color the
proteins based on their protein IDs. If the protein ID matches the
protein ID of the model, then just give the the color specified in
the model. This requires determining if the species in the analysis
is the same species as the model.

The second way to color the proteins is to use the HMMs of them to
find the best hits, then color with the specified color in the rbh
file accompanying the hmms.
"""

import os
import pandas as pd
import sys

#def find_matching_species(spA_prot_list, LG_rbh_path):
#    """
#    """

class LG_db:
    """
    This class ingests a directory of an odp linkage database,
     checks to see if the input is legal based on the LG db spec,
     and makes data structures to access the information later.
    """
    def __init__(self, LG_db_name, LG_db_directory,
                       hmm_results_paths):
        self.name      = LG_db_name
        # check if directory is OK.
        if not os.path.isdir(LG_db_directory):
            raise IOError("{} for LG database {} does not exist".format(LG_db_directory, LG_db_name))
        self.directory = os.path.abspath(LG_db_directory)
        self._dirfiles = list(os.listdir(self.directory))

        # check that there is only one .rbh file in the directory
        if len([x for x in self._dirfiles if x.endswith(".rbh")]) != 1:
            raise IOError("There must be a single .rbh file in the LG db directory for {}.\n Instead we found {}".format(
                LG_db_name, [x for x in self.dirfiles if x.endswith(".rbh")]))
        self.rbhfile = os.path.join(self.directory, [x for x in self._dirfiles if x.endswith(".rbh")][0])

        # check that there is only one .hmm file in the directory
        if len([x for x in self._dirfiles if x.endswith(".hmm")]) != 1:
            raise IOError("There must be a single .hmm file in the LG db directory for {}.\n Instead we found {}".format(
                LG_db_name, [x for x in self.dirfiles if x.endswith(".hmm")]))
        self.hmmfile = os.path.join(self.directory, [x for x in self._dirfiles if x.endswith(".hmm")][0])

        # now open the rbh file and build appropriate data structures
        self.rbhdf = pd.read_csv(self.rbhfile, sep = "\t")
        self.rbhspecies = [x.strip("_scaf") for x in self.rbhdf.columns if x.endswith("_scaf")]

        # check that the requisite columns are in the rbh file
        for thiscol in ["rbh", "gene_group", "color"]:
            if thiscol not in self.rbhdf:
                raise IOError("The rbh df must have the following column, but does not currently: {}".format(
                    thiscol))

        # make tables that are easy to look up gene and species to color and group
        self.sp_to_gene_to_color = self._gen_sp_to_gene_to_color()
        self.sp_to_gene_to_group = self._gen_sp_to_gene_to_group()
        self.rbh_to_color = dict(zip(self.rbhdf["rbh"], self.rbhdf["color"]))
        self.rbh_to_group = dict(zip(self.rbhdf["rbh"], self.rbhdf["gene_group"]))

        # Open the hmm file and parse it to check that each ortholog has a HMM.
        rbh_orthologs = set(list(self.rbhdf["rbh"]))
        hmm_orthologs = set()
        with open(self.hmmfile, "r") as f:
            for line in f:
                line = line.strip()
                if line and line.startswith("NAME "):
                    hmm_orthologs.add(line.replace("NAME ", "").strip())

        # if any of these fail, the rbh and hmm files don't match completely
        if len(rbh_orthologs - hmm_orthologs) != 0:
            raise IOError("The rbh file contains some orthologs not in the hmm file: {}".format(
                rbh_orthologs - hmm_orthologs))
        if len(hmm_orthologs - rbh_orthologs) != 0:
            raise IOError("The hmm file contains some orthologs not in the rbh file: {}".format(
                rbh_orthologs - hmm_orthologs))

        # now make a dict of the most likely colors for each group
        self.group_to_color_df = self.rbhdf.groupby(
                 by=["gene_group"])["color"].value_counts(
                 ascending = False).rename("Counts").reset_index(
                 ).drop_duplicates("gene_group")
        self.group_to_color = dict(zip(self.group_to_color_df["gene_group"],
                                        self.group_to_color_df["color"]))
        for this in ["None", "none", "na", "NA"]:
            self.group_to_color[this] = "#000000"

        # parse the hmm results. make a dict of proteins to the rbh
        self.hmm_prot_to_rbh = self._gen_hmm_prot_to_rbh(hmm_results_paths)

        # keep track of how we assigned colors, this is useful later
        self.color_method = ""

    def _gen_hmm_prot_to_rbh(self, hmm_results_paths):
        """
        Read in the HMM results. They must be in blastp format.
        column 1 is the rbh that we will look for.
        column 2 is the protein in the species' genomes being plotted
        second to last column is evalue.
        last column is bitscore
        """
        dataframes = []
        for thisfile in hmm_results_paths:
            if not os.path.exists(thisfile):
                raise IOError("{} does not exist.".format(thisfile))
            dataframes.append(pd.read_csv(thisfile, sep = "\t", header = None))
            dataframes[-1].columns =  ["qseqid", "sseqid", "pident", "length",
                                       "mismatch", "gapopen", "qstart", "qend",
                                       "sstart", "send", "evalue", "bitscore"]
        done = False;
        keeps = []
        df = pd.concat(dataframes).sort_values(
              by=["evalue", "bitscore"], ascending = [True, False]
              ).reset_index(drop=True)
        df = df.loc[df["evalue"] <= 1e-5,]

        # get the best hit for each
        while len(df) > 0:
            qseqid = df.iloc[0]["qseqid"]
            sseqid = df.iloc[0]["sseqid"]
            keeps.append(df.iloc[0].to_dict())
            df = df.loc[df["qseqid"] != qseqid,]
            df = df.loc[df["sseqid"] != sseqid,]
        df = pd.DataFrame.from_dict(keeps)
        if not len(df["qseqid"].unique()) == len(df):
            raise IOError("We should not see this, qseqid")
        if not len(df["sseqid"].unique()) == len(df):
            raise IOError("We should not see this, sseqid")
        return dict(zip(df["sseqid"],
                        df["qseqid"]))

    def _gen_sp_to_gene_to_color(self):
        """
        Return a dictionary that makes it easier to look up what color corresponds
         to what gene.
        """
        tempdict = {sp: {} for sp in self.rbhspecies}
        for index, row in self.rbhdf.iterrows():
            for thissp in self.rbhspecies:
                color = row["color"] if "#" in row["color"] else "#000000"
                tempdict[thissp][row["{}_gene".format(thissp)]] = color
        return tempdict

    def _gen_sp_to_gene_to_group(self):
        """
        Return a dictionary that makes it easier to look up what group corresponds
         to what gene.
        """
        tempdict = {sp: {} for sp in self.rbhspecies}
        for index, row in self.rbhdf.iterrows():
            for thissp in self.rbhspecies:
                tempdict[thissp][row["{}_gene".format(thissp)]] = row["gene_group"]
        return tempdict

    def _sp_matches_which_db_species(self, protein_id_list):
        """
        As input, takes a list of protein ids. This method uses this list
         of protein IDs to determine if the input protein IDs correspond to
         a particular species in the LG database.

        The method for doing this is checking that there is at least a 10:1
         ratio of the proteins in the protein_id_list corresponding to
         one species over another in the .rbh file.
        """
        species_count = {key: 0 for key in self.rbhspecies}
        for thissp in self.rbhspecies:
            species_count[thissp] =  self.rbhdf["{}_gene".format(thissp)].isin(protein_id_list).sum()

        species_count = {key: species_count[key] for key in species_count
                         if species_count[key] != 0}

        # if there are no matches, then just return none
        if len(species_count) == 0:
            return None
        else:
            # there may be a match
            if len(species_count) == 1:
                return [k for k in species_count][0]
            else:
                # get the largest
                species_count_sorted_keys = [[k,v] for k, v in sorted(
                    species_count, key=species_count.get, reverse=True).items()]
                # if the first is at least 10 times larger than the second, return it
                if species_count_sorted_keys[0][1] >= (species_count_sorted_keys[0][2] * 10):
                    return species_count_sorted_keys[0][0]
                else:
                    return None
        # catch case
        return None

    def color_dataframe(self, plotdf):
        """
        - The plotdf is the dataframe that we wish to color
        """
        # First step, determine if the species in the plotdf are the same
        #  as in the LG database
        all_species = [x.strip("_scaf") for x in plotdf.columns if x.endswith("_scaf")]

        # check which species in plotdf correspond with which species in
        species_to_LG_species = {}
        for thissp in all_species:
            temp = self._sp_matches_which_db_species(plotdf["{}_gene".format(thissp)])
            if temp: # if not None
                species_to_LG_species[thissp] = temp

        if len(species_to_LG_species) > 0:
            self.color_method = "protein ids of {} in the {} models".format(
                "+".join([x for x in species_to_LG_species]), self.name)

            for index, row in plotdf.iterrows():
                color = "#000000"
                group = "None"
                for thissp in species_to_LG_species:
                    plotgene = row["{}_gene".format(thissp)]
                    othersp = species_to_LG_species[thissp]
                    if plotgene in self.sp_to_gene_to_color[othersp]:
                        thiscolor = self.sp_to_gene_to_color[othersp][plotgene]
                        if ("#" in thiscolor) and (thiscolor != "#000000"):
                            color = thiscolor
                    if plotgene in self.sp_to_gene_to_group[othersp]:
                        thisgroup = self.sp_to_gene_to_group[othersp][plotgene]
                        if thisgroup != "None":
                            group = thisgroup
                plotdf.loc[index,"color"] = color
                plotdf.loc[index,"gene_group"] = group
        else:
            # There is no species match, so we should use HMMs to find the proteins
            self.color_method = "HMM of {} models against the sp-sp proteins".format(self.name)
            for index, row in plotdf.iterrows():
                color = "#000000"
                group = "None"
                for thissp in all_species:
                    plotgene = row["{}_gene".format(thissp)]
                    if plotgene in self.hmm_prot_to_rbh:
                        this_rbh = self.hmm_prot_to_rbh[plotgene]
                        thisgroup = self.rbh_to_group[this_rbh]
                        thiscolor = self.rbh_to_color[this_rbh]
                        if thisgroup != "None":
                            group = thisgroup
                        if ("#" in thiscolor) and (thiscolor != "#000000"):
                            color = thiscolor
                plotdf.loc[index,"color"] = color
                plotdf.loc[index,"gene_group"] = group

        return plotdf
