#!/usr/local/bin python3

from ete3 import Tree, TreeStyle, AttrFace, faces, NodeStyle, NCBITaxa

import sys

def main():
    # Read in the Newick tree from the first sys arg.
    # Parse the node names and add them to the tree. These are NCBI TaxIDs.
    tree = Tree(sys.argv[1], format=1)

    # This is the frequency of the different ranks in a dataset of about 3000 species
    # clade 16144
    # no rank 3312
    # superkingdom 1446
    # kingdom 1446
    # phylum 1438
    # class 1391
    # order 1285
    # subphylum 1252
    # family 996
    # subclass 970
    # infraclass 931
    # suborder 900
    # cohort 776
    # superfamily 701
    # infraorder 621
    # superclass 614
    # subfamily 594
    # superorder 550
    # genus 483
    # parvorder 214
    # tribe 212
    # subgenus 67
    # subcohort 58
    # species 37
    # species group 22
    # species subgroup 12
    # subtribe 9
    # section 1
    # series 1
    # start the taxonomy lookup
    ncbi = NCBITaxa()
    # Collapse all of the nodes that are "Class" level in NCBI taxonomy.
    #  The node IDs are NCBI TaxID, so we can use that to identify the nodes.
    node_name_transform = {}
    collapse_list = set()
    for node in tree.traverse():
        if node.name.isdigit() and int(node.name) == 6605:
            print("We found cephalopods!!")
            sys.exit()
        # Check that the node name isn't empty, and that the node name is a parsable integer.
        # If the node name isn't a parseable integer
        if node.name != "" and node.name.isdigit():
            # lookup node name in NCBI taxonomy. The node name
            NCBIname = ncbi.get_taxid_translator([node.name])
            # Modify the node_name_transform dictionary to include the NCBI node name.
            # We want the new node name to be the NCBI taxid, underscore, and the NCBI name.
            node_name_transform[node.name] = "{}_{}".format(node.name, NCBIname[int(node.name)])
            lineage = ncbi.get_lineage(node.name)
            rank = ncbi.get_rank(lineage)
            if "ephal" in NCBIname[int(node.name)]:
                print(NCBIname[int(node.name)], lineage)
            if rank[lineage[-1]] == "class":
                collapse_list.add(node.name)
    # Collapse the nodes in the collapse list.
    tree.prune(collapse_list, preserve_branch_length=True)

    # For the remaining nodes, modify the name using the node_name_transform dictionary.
    for node in tree.traverse():
        if node.name in node_name_transform:
            node.name = node_name_transform[node.name]


    # Define a custom style for nodes
    nstyle = NodeStyle()
    nstyle["size"] = 0

    # Apply the style to collapsed nodes
    for node in tree.traverse():
        if not node.is_leaf():
            node.set_style(nstyle)

    # Create a TreeStyle with a specific layout
    ts = TreeStyle()
    ts.show_leaf_name = True
    ts.show_branch_length = False  # Hide branch lengths for clarity

    # Add additional information along each branch
    def my_layout(node):
        if node.is_leaf():
            name_face = AttrFace("name", fsize=10)
            faces.add_face_to_node(name_face, node, column=0)

    # Apply the layout to the tree
    ts.layout_fn = my_layout

    # Render the tree to a file (e.g., "cladogram.pdf")
    tree.render("cladogram.pdf", tree_style=ts)

if __name__ == '__main__':
    main()