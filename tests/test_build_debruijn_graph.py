#!/usr/bin/env python
import unittest
import sys
sys.path.append('../src')
from contigs import build_debruijn_graph

class TestDBG(unittest.TestCase):

    def test_build_debruijn_graph(self):

        id_presuffix_dict = {'>test_seq1': [['AT', 'TC'],['TC', 'CA'],['CA', 'AG'],['AG', 'GG'],['GG', 'GT'],['GT', 'TA'],['TA', 'AC']],
                                '>test_seq2': [['GT', 'TA'], ['TA', 'AA'], ['AA', 'AC']]}
        exp_edges_hash={'AT': ['TC'],'TC': ['CA'],
                        'CA': ['AG'],'AG': ['GG'],
                        'GG': ['GT'],'GT': ['TA', 'TA'],
                        'TA': ['AC', 'AA'],'AA': ['AC']}

        edges_hash=build_debruijn_graph(id_presuffix_dict)

        self.assertEqual(edges_hash,exp_edges_hash)

if __name__ == '__main__':
    unittest.main()
