#!/usr/bin/env python
import unittest
import sys
sys.path.append('../src')
from contigs import calculate_node_degree

class TestNodeDegree(unittest.TestCase):

    def test_calculate_node_degree(self):

        id_presuffix_dict = {'>test_seq1': [['AT', 'TA'], ['TA', 'AG'],['AG','GT'],['GT','TC'],['TC','CA'],['CA','AT']], 
                                '>test_seq2': [['TT', 'TG'],['TG','GC'],['GC','CT'],['CT','TA'],['TA','AT']]}
        
        exp_node_degree = {'AT': [2, 1],'TA': [2, 2],
            'AG': [1, 1],'GT': [1, 1],
            'TC': [1, 1],'CA': [1, 1],
            'TT': [0, 1],'TG': [1, 1],
            'GC': [1, 1],'CT': [1, 1]}
        
        exp_start_nodes=['TT']
        exp_end_nodes=['AT']

        nodes_degree, starting_nodes, ending_nodes = calculate_node_degree(id_presuffix_dict)

        self.assertEqual(nodes_degree, exp_node_degree)
        self.assertEqual(starting_nodes, exp_start_nodes)
        self.assertEqual(ending_nodes, exp_end_nodes)

if __name__ == '__main__':
    unittest.main()    