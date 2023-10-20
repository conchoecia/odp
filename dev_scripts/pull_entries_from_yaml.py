#!/usr/bin/env python

# This takes a list of entries that are expected to be in a yaml file, plus a yaml file,
# then prints out only the desired entries

import argparse
import yaml

def parse_args():
    """
    We need two args:
      - keeps: a list of entries to keep that we expect to appear in the yaml file
      - yaml:  the yaml file from which we will pull entries. This is formatted like the other yaml files for odp
    """
    parser = argparse.ArgumentParser(description='Pull entries from yaml file')
    parser.add_argument('-k', '--keeps', nargs='+', help='list of entries to keep')
    parser.add_argument('-y', '--yaml', help='yaml file to pull entries from')
    args = parser.parse_args()
    return args

def get_genus_species_list(keeps_file):
    """
    Gets the species combos to keep from the keeps file.
    Each line in the keeps file is the species string of a single file.
    Lines that begin with # characters are ignored.
    Returns a list of species strings.
    """
    keeps_species = []
    for thisfile in keeps_file:
        with open(thisfile, 'r') as f:
            for line in f:
                line = line.strip()
                if line:
                    if line.startswith('#'):
                        continue
                    keeps_species.append(line.strip())
    return keeps_species

def get_yaml_entries(yaml_file):
    """
    This opens the yaml file and loads it as a dictionary
    """
    yaml_dict = {}
    with open(yaml_file, 'r') as f:
        #I was getting this error TypeError: load() missing 1 required positional argument: 'Loader'
        # so I can't use         yaml_dict = yaml.load(f)
        yaml_dict = yaml.load(f, Loader=yaml.FullLoader)
    return yaml_dict

def filter_yaml_entries(yaml_dict, keeps_species):
    """
    Just keeps the entries that are matched in keeps_species,

    We can either match on exact matches to the sample name or an exact match to the species name
    """
    keep_these_entries = {}

    # we will use this dict for genus species matches later on.
    genus_species = {}
    for thisstring in keeps_species:
        fields = thisstring.split(" ")
        if len(fields) > 1:
            # we have found something that could be a genus species combo
            thisgenus   = fields[0]
            thisspecies = " ".join(fields[1::])
            if thisgenus not in genus_species:
                genus_species[thisgenus] = set()
            genus_species[thisgenus].add(thisspecies)
    genus_species = {k: list(v) for k, v in genus_species.items()}

    # check if species is in the yaml dict
    if not "species" in yaml_dict:
        raise IOError("The species field wasn't found in the yaml file, but it should have been.")

    # now we go through the dict and pick out the species we want to keep
    for thissp in yaml_dict["species"]:
        # first we check to see if the species is in the keep list.
        # If so, this means the user put the literal species sample name in the keep list and we don't have to do anything else
        if thissp in keeps_species:
            keep_these_entries[thissp] = yaml_dict["species"][thissp]
        else:
            # there wasn't an exact match, so now we have to do some more advanced matching
            # Check if there was a genus and species defined in this sample's entry in the yaml file
            if "genus" in yaml_dict["species"][thissp] and "species" in yaml_dict["species"][thissp]:
                # now we check if the genus is in the keep list
                if yaml_dict["species"][thissp]["genus"] in genus_species:
                    # now we check if the species is in the keep list
                    if yaml_dict["species"][thissp]["species"] in genus_species[yaml_dict["species"][thissp]["genus"]]:
                        keep_these_entries[thissp] = yaml_dict["species"][thissp]

    return keep_these_entries

def main():
    args = parse_args()
    keep_species = get_genus_species_list(args.keeps)
    yaml_dict = get_yaml_entries(args.yaml)

    # now we filter out entries that are not in the keep list
    keep_these = filter_yaml_entries(yaml_dict, keep_species)

    # now we print out a text file of these entries. Dump as a yaml to std.out
    species_dict = {"species": keep_these}
    print(yaml.dump(species_dict, default_flow_style=False))

if  __name__ == "__main__":
    main()