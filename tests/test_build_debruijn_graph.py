#!/usr/bin/env python
import unittest
import sys
sys.path.append('../src')
from contigs import build_debruijn_graph

class TestDBG(unittest.TestCase):

    def test_build_debruijn_graph(self):

        id_presuffix_dict = {'>test_seq1': [['AT', 'TA'], ['TA', 'AG'],['AG','GT'],['GT','TC'],['TC','CA'],['CA','AT']], 
                                '>test_seq2': [['TT', 'TG'],['TG','GC'],['GC','CT'],['CT','TA'],['TA','AT']]}
        
        exp_edges_hash={
            'AT': ['TA'],'TA': ['AG', 'AT'],
            'AG': ['GT'], 'GT': ['TC'],
            'TC': ['CA'],'CA': ['AT'],
            'TT': ['TG'],'TG': ['GC'],
            'GC': ['CT'],'CT': ['TA']}

        edges_hash=build_debruijn_graph(id_presuffix_dict)

        self.assertEqual(edges_hash,exp_edges_hash)

if __name__ == '__main__':
    unittest.main()
