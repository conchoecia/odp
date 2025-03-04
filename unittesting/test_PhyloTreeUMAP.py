# insert the path to the source file
# ODP-specific imports
import sys
import os

thisdir = os.path.dirname(os.path.abspath(__file__))
source_path = os.path.join(thisdir, "../dev_scripts")
if source_path not in sys.path:
    sys.path.insert(1, source_path)
import PhyloTreeUMAP


import unittest   # The test framework


class Test_taxids_to_analyses(unittest.TestCase):
    def test_empty(self):
        """
        No analyses are passed. Should fail with an error.
        """
        analysis = []
        self.assertRaises(ValueError, PhyloTreeUMAP.taxids_to_analyses, analysis)

    def test_entry_empty(self):
        """
        In this case, there is an entry that has nothing in it. This isn't allowed, as we require some taxid to look at.
        For example - the first one is fine, the second should trigger an error.
          [
           [[6340],  [42113]],
           [[],  [42113]]
          ]
        """
        analysis = [ [[6340],  [42113]],
                     [[],  [42113]]
                   ]
        self.assertRaises(ValueError, PhyloTreeUMAP.taxids_to_analyses, analysis)

    def test_entry_empty2(self):
        """
        Makes sure that each entry is length two. The two things are two lists. [[6340],  [42113]]
        """
        analysis = [[6340],  [42113], []]
        # should fail because the length is 3
        self.assertRaises(ValueError, PhyloTreeUMAP.taxids_to_analyses, analysis)
        # should fail because one is not a list
        analysis = [[6340],  42113]
        self.assertRaises(ValueError, PhyloTreeUMAP.taxids_to_analyses, analysis)
        # should fail because another is not a list
        analysis = [{},  [42113]]
        self.assertRaises(ValueError, PhyloTreeUMAP.taxids_to_analyses, analysis)

    def test_integers(self):
        """
        Make sure that the taxids are all integers. This will make it fail.
        """
        analysis = [[6340, "ShouldFail"],  []]
        self.assertRaises(ValueError, PhyloTreeUMAP.taxids_to_analyses, analysis)

    def test_annelida(self):
        """
        Just the annelida
        """
        analysis = [ [[6340],  []] ]
        results = {"Annelida_6340_without_None": [[6340],  []]}
        self.assertEqual(PhyloTreeUMAP.taxids_to_analyses(analysis), results)

    def test_multi_without_Beroida(self):
        """
        Should have multiple clades included, one excluded
        """
        analysis = [ [[10197, 6040, 6073], [37538]]] #ctenos, sponges, cnidarians
        results = {"Porifera_Cnidaria_Ctenophora_6040_6073_10197_without_37538": [[10197, 6040, 6073], [37538]]}
        self.assertEqual(PhyloTreeUMAP.taxids_to_analyses(analysis), results)
        # should be the same thing if we shuffle the input
        analysis = [ [[6073, 10197, 6040], [37538]]] #ctenos, sponges, cnidarians
        results = {"Porifera_Cnidaria_Ctenophora_6040_6073_10197_without_37538": [[6073, 10197, 6040], [37538]]}
        self.assertEqual(PhyloTreeUMAP.taxids_to_analyses(analysis), results)


if __name__ == '__main__':
    unittest.main()