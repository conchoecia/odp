"""
Functions related to Fisher's Exact Test (FET).

Does not appear to need additional dependencies.
"""

def calculate_FET(genesdf, scafdict):
    """
    Takes in a rbh file and calculates the false discovery rate using
      Fisher's Exact Test.

    Scaffolds that are not included in scafdict are removed from the analysis

    Inputs:
      - genesdf - a pandas dataframe of RBH
      - species_list - a list of species?? why is this here
      - minscafsizelist

    Failure cases:
      - there are more than two species in the genesdf
    """
    import scipy.stats as stats
    all_species = [x.replace("_gene", "") for x in genesdf.columns if x.endswith("_gene")]
    # make sure there are only two species since this method only works with that
    if len(all_species) > 2:
        raise IOError("This method can only be performed with two species. This has {}".format(all_species))

    # now filter rows based on the scaffolds that we allow
    for i in [0,1]:
        # make sure this species is in the scafdict keys
        if all_species[i] not in scafdict:
            raise IOError("This species is not in the scafdict: {}".format(all_species[i]))
        # we use a set because this is better for filtering
        allowed_scaffolds = set(scafdict[all_species[i]])  # Convert to a set for fast lookup
        scafcol = "{}_scaf".format(all_species[i])
        genesdf = genesdf.loc[genesdf[scafcol].isin(allowed_scaffolds)]
    #print("genesdf after filtering\n", genesdf)

    genesdf["whole_FET"] = -1
    genesdf["break_FET"] = -1
    # now perform Fisher's Exact Test
    if len(all_species) == 2:
        sp1 = list(sorted(all_species))[0]
        sp2 = list(sorted(all_species))[1]
        for sp1scaf, sp2scaf, FET in [
                ["{}_scaf".format(sp1), "{}_scaf".format(sp2), "whole_FET"],
                ["{}_breakchrom".format(sp1), "{}_breakchrom".format(sp2), "break_FET"]]:
            # make this table for the chroms
            sp1_chroms = genesdf[sp1scaf].unique()
            sp2_chroms = genesdf[sp2scaf].unique()
            combos = {(x, y): 1 for x in sp1_chroms for y in sp2_chroms}
            for thiscombo in combos:
                sp1_chrom = thiscombo[0]
                sp2_chrom = thiscombo[1]
                in1in2  =  len(genesdf.loc[(genesdf[sp1scaf] == sp1_chrom) & (genesdf[sp2scaf] == sp2_chrom),])
                out1in2 =  len(genesdf.loc[(genesdf[sp1scaf] != sp1_chrom) & (genesdf[sp2scaf] == sp2_chrom),])
                in1out2 =  len(genesdf.loc[(genesdf[sp1scaf] == sp1_chrom) & (genesdf[sp2scaf] != sp2_chrom),])
                out1out2 = len(genesdf.loc[(genesdf[sp1scaf] != sp1_chrom) & (genesdf[sp2scaf] != sp2_chrom),])
                table = [[in1in2, out1in2], [in1out2, out1out2]]
                oddsratio, pvalue = stats.fisher_exact(table, alternative="greater")
                combos[thiscombo] = pvalue * len(combos)
            for index, row in genesdf.iterrows():
                rowsp1 = row[sp1scaf]
                rowsp2 = row[sp2scaf]
                genesdf.loc[index, FET] = combos[(row[sp1scaf], row[sp2scaf])]
        #genesdf = genesdf.groupby(["{}_scaf".format(sp1), "{}_scaf".format(sp2)])

    return genesdf


def FET_bfs(graph, start):
    """
    This code taken from this SO question: https://stackoverflow.com/questions/53573865
     and this specific answer: https://stackoverflow.com/a/53574094/5843327
     by SO user: https://stackoverflow.com/users/4001592/dani-mesejo
    """
    visited, queue = set(), [start]
    while queue:
        vertex = queue.pop(0)
        if vertex not in visited:
            visited.add(vertex)
            queue.extend(graph[vertex] - visited)
    return visited

def FET_connected_components(G):
    """
    This code taken from this SO question: https://stackoverflow.com/questions/53573865
     and this specific answer: https://stackoverflow.com/a/53574094/5843327
     by SO user: https://stackoverflow.com/users/4001592/dani-mesejo
    """
    seen = set()
    for v in G:
        if v not in seen:
            c = set(FET_bfs(G, v))
            yield c
            seen.update(c)

def FET_graph(edge_list, G = None):
    """
    This code taken from this SO question: https://stackoverflow.com/questions/53573865
     and this specific answer: https://stackoverflow.com/a/53574094/5843327
     by SO user: https://stackoverflow.com/users/4001592/dani-mesejo

    The parameters are set up so that a graph can be made from just a list of edges,
        or a graph can be added to an existing graph.
    """
    result = {}
    if G is not None:
        result = G
    for source, target in edge_list:
        result.setdefault(source, set()).add(target)
        result.setdefault(target, set()).add(source)
    return result

def FET_sample_to_chrom_to_CC(components, xsample, ysample, splitter):
    """
    Turns this into something we can use to assign connected components.
    """
    # convert components into a dictionary for the x or y-chroms lookup. We just want the CC, because we can filter on these later
    sample_to_chrom_to_CC = {xsample: {}, ysample: {}}
    cc_counter = 0
    for CC in components:
        for scaf in CC:
            esample = scaf.split(splitter)[0]
            echrom  = scaf.split(splitter)[1]
            # check if we have seen this already, if so, raise an error
            if sample_to_chrom_to_CC[esample].get(echrom, None) is not None:
                raise ValueError("ERROR: we have seen this chrom before, but it is in a different CC. This should not happen. We saw {}.".format(echrom))
            sample_to_chrom_to_CC[esample][echrom] = cc_counter
        cc_counter += 1
    return sample_to_chrom_to_CC