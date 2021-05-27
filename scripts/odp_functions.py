"""
These are functions that are shared by odp, odp_trio, and other scripts
"""
from itertools import combinations
from itertools import product

def check_legality(config):# check for legal config entries. Useful fr finding misspelled entries
    legal = ["proteins", "prot_to_loc", "genome", "minscafsize", "manual_breaks", "chrom_to_color", "plotorder"]
    illegal = set()
    for this_axis in ["xaxisspecies", "yaxisspecies"]:
        if this_axis in config:
            for this_sample in config[this_axis]:
                for key in config[this_axis][this_sample]:
                    if key not in legal:
                        illegal.add(key)
    if len(illegal) > 0:
        print("We found some fields in your config file that are not used by this program.")
        print("The only fields allowed for individual samples are:")
        for key in legal:
            print("  - {}".format(key))
        print("The keys that we found that are not allowed/in the list above are:")
        for key in illegal:
            print("  - {}".format(key))
        sys.exit()

def flatten(list_of_lists):
    """flatten a list of lists, unique only"""
    return list(set([item for sublist in list_of_lists for item in sublist]))

def expand_avoid_matching_x_and_y(filestring, xsamples, ysamples):
    """
    this works like snakemake's expand function but does not generate
     files where xsample equals ysample

    outputs a list of files
    """
    outlist = []
    for xsamp in xsamples:
        for ysamp in ysamples:
            if xsamp != ysamp:
                outlist.append(filestring.format(xsamp, ysamp))
    return outlist

def expand_avoid_matching_x_and_y_third(filestring, xsamples, ysamples, third):
    """
    this works like snakemake's expand function but does not generate
     files where xsample equals ysample

    outputs a list of files.

    Note March 20th, 2021 - the "third" aspect is related to coloring by another species
    """
    outlist = []
    for xsamp in xsamples:
        for ysamp in ysamples:
            if xsamp != ysamp:
                for t in third:
                    print("hey", xsamp, ysamp, third)
                    outlist.append(filestring.format(xsamp, ysamp, t))
    return outlist

def expand_avoid_matching(filestring, **kwargs):
    """
    this works like snakemake's expand function but does not generate
     files where xsample equals ysample

    outputs a list of files
    """
    outlist = []
    keys = [x for x in kwargs]
    values = [[y for y in kwargs[x]] for x in keys]
    nonmatching_products = [x for x in list(product(*values)) if len(set(x)) == len(x)]
    for entry in nonmatching_products:
        these_kwargs = {}
        for i in range(len(entry)):
            these_kwargs[keys[i]] = entry[i]
        outlist.append(filestring.format(**these_kwargs))
    return [x for x in outlist]

def generate_coord_structs_from_chrom_to_loc(prot_to_loc_file):
    """
    This parses a .chrom file and outputs five data structures that are easily
     used for mapping pandas dataframes.
    The output is a dict of dicts. Not the most intuitive format but easy for
     mapping to column values.
     { "prot_to_scaf":   prot_to_scaf,
       "prot_to_strand": prot_to_strand,
       "prot_to_start":  prot_to_start,
       "prot_to_stop":   prot_to_stop,
       "prot_to_middle": prot_to_middle }
    """
    prot_to_scaf   = {}
    prot_to_strand = {}
    prot_to_start  = {}
    prot_to_stop   = {}
    prot_to_middle = {}
    print("prot_to_loc_file", prot_to_loc_file)
    with open(prot_to_loc_file, "r") as f:
       for line in f:
           line = line.strip()
           if line:
               splitd = line.split()
               prot = splitd[0]
               # add things now
               prot_to_scaf[prot]   = splitd[1]
               prot_to_strand[prot] = splitd[2]
               start = int(splitd[3])
               prot_to_start[prot]  = start
               stop = int(splitd[4])
               prot_to_stop[prot]   = stop
               stop = int(splitd[4])
               prot_to_middle[prot] = int(start + (stop - start)/2)
    return { "prot_to_scaf":   prot_to_scaf,
             "prot_to_strand": prot_to_strand,
             "prot_to_start":  prot_to_start,
             "prot_to_stop":   prot_to_stop,
             "prot_to_middle": prot_to_middle }
