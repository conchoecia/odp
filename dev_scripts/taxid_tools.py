"""
These are functions that are used to work with ncbi taxids.
"""

import re
import warnings

def NCBI_taxid_to_taxdict(ncbi, taxid) -> dict:
    """
    Takes a single NCBI taxid as input and returns a dictionary with useful information:

    Input:
      - ncbi:  The NCBITaxa object
      - taxid: The NCBI taxid
    Output:
      - A dictionary with the following
        taxid: The taxid, same as the input
        taxname: The name of this specific taxid
        taxname_list: A list of the taxonomy names, like ["root", "cellular organisms", "Eukaryota", "Opisthokonta"]
        taxid_list: A list of the taxids, like [1, 131567, 2759, 33154]
        level_1: The first level of the taxid, like "root (1); cellular organisms (131567); Eukaryota (2759); Opisthokonta (33154)"
        ... up to level_10
        printstring: The printstring of the taxid, like "root (1); cellular organisms (131567); Eukaryota (2759); Opisthokonta (33154)"
    """
    if isinstance(taxid, str):
        # first check that the taxid is an integer
        if not re.match(r"^[0-9]*$", taxid):
            raise ValueError(f"There is a non-numeric character in the taxid string, {taxid}, for file {thisfile}. Exiting.")
    elif isinstance(taxid, int):
        pass
    else:
        raise ValueError(f"The taxid is not a string or an integer. It is a {type(taxid)}. Exiting.")

    original_taxid = taxid
    # now we fix the taxid if it is something that existed, but no longer existed.
    old_translator = {876063: 3126489, # this is for the moth Ochlodes sylvanus
                      355208: 3056719, # this is for the moth Spicauda simplicius
                      }
    if taxid in old_translator:
        taxid = old_translator[taxid]

    # safe, get the lineage
    entry = {"taxid": taxid}
    # for each node, make the full lineage string, in this form "Metazoa;Bilateria;Protostomes"
    # If there is a taxid change we need to catch the warning warnings.warn("taxid %s was translated into %s" %(taxid, merged_conversion[taxid]))
    #  We don't care if they changed the taxid name from what was originally in the file.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        lineage = ncbi.get_lineage(taxid)
    # We're just going to ignore whatever changes were made, and we will force the last entry of the lineage
    #  to be what is in the filename.
    names = ncbi.get_taxid_translator(lineage)
    names[original_taxid] = names[lineage[-1]]
    lineage[-1] = original_taxid
    # ^^ now we're done fidgeting with the lineage and it should work even if the lineage has been changed since the genome was downloaded.
    # make sure that the lineage and the names are not empty
    if len(lineage) == 0:
        raise ValueError(f"The lineage is empty for the taxid {taxid}. Exiting.")
    if len(names) == 0:
        raise ValueError(f"The names are empty for the taxid {taxid}. Exiting.")
    entry["taxname"]          = names[taxid]
    entry["taxid_list"]       = [taxid for taxid in lineage]
    entry["taxid_list_str"]   = ";".join([str(taxid) for taxid in lineage])
    entry["taxname_list"]     = [names[taxid] for taxid in lineage]
    entry["taxname_list_str"] = ";".join([names[taxid] for taxid in lineage])

    npl = 4 # number of taxonomic units per level
    num_rows = 10 # what level do we want to go to
    for i in range(1, num_rows+1):
        thislevel = f"level_{i}"
        j = (i-1)*npl
        entry[thislevel] = ";".join([f" {names[taxid]} ({taxid})" for taxid in lineage[j:j+npl]])

    entry["printstring"]  = ";".join([f" {names[taxid]} ({taxid})" for taxid in lineage])
    return entry
