#!/usr/bin/env python
import unittest
import sys
sys.path.append('../src')
from contigs import calculate_node_degree

class TestNodeDegree(unittest.TestCase):

    def test_calculate_node_degree(self):

        id_presuffix_dict = {'>test_seq1': [['AT', 'TC'],['TC', 'CA'],['CA', 'AG'],['AG', 'GG'],['GG', 'GT'],['GT', 'TA'],['TA', 'AC']],
                                '>test_seq2': [['GT', 'TA'], ['TA', 'AA'], ['AA', 'AC']]}
        
        exp_node_degree = {'AT': [0, 1],
                            'TC': [1, 1],
                            'CA': [1, 1],
                            'AG': [1, 1],
                            'GG': [1, 1],
                            'GT': [1, 2],
                            'TA': [2, 2],
                            'AC': [2, 0],
                            'AA': [1, 1]}
        
        exp_start_nodes=['AT','GT']
        exp_end_nodes=[]

        nodes_degree, starting_nodes, ending_nodes = calculate_node_degree(id_presuffix_dict)

        self.assertEqual(nodes_degree, exp_node_degree)
        self.assertEqual(starting_nodes, exp_start_nodes)
        self.assertEqual(ending_nodes, exp_end_nodes)

if __name__ == '__main__':
    unittest.main()    