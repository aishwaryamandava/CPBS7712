#!/usr/bin/env python
import unittest
import sys
sys.path.append('../src')
from contigs import generate_kmer_nodes 

class TestKmer(unittest.TestCase):

    def test_generate_kmer_nodes(self):

       seq_path='./test_READS.fasta'
       k_len=3

       id_kmer_dict, id_read_dict, id_presuffix_dict, nodes_unique, edges_all, suffixes, kmer_id_coord = generate_kmer_nodes(seq_path, k_len)
    
       exp_id_kmer_dict = {'>test_seq1': ['ATC', 'TCA', 'CAG', 'AGG', 'GGT', 'GTA', 'TAC'],
                           '>test_seq2': ['GTA', 'TAA', 'AAC']}
       exp_id_read_dict = {'>test_seq1': 'ATCAGGTAC', '>test_seq2': 'GTAAC'}
       exp_id_presuffix_dict = {'>test_seq1': [['AT', 'TC'],['TC', 'CA'],['CA', 'AG'],['AG', 'GG'],['GG', 'GT'],['GT', 'TA'],['TA', 'AC']],
                                '>test_seq2': [['GT', 'TA'], ['TA', 'AA'], ['AA', 'AC']]}
       exp_nodes_unique = {'AA', 'AC', 'AG', 'AT', 'CA', 'GG', 'GT', 'TA', 'TC'}
       exp_edges_all = [('AT', 'TC'),('TC', 'CA'),('CA', 'AG'),('AG', 'GG'),('GG', 'GT'),('GT', 'TA'),('TA', 'AC'),('GT', 'TA'),('TA', 'AA'),('AA', 'AC')]
       exp_suffixes = {'AA', 'AC', 'AG', 'CA', 'GG', 'GT', 'TA', 'TC'}
       exp_kmer_id_coord = {'ATC': {'>test_seq1': [(0, 2)]},'TCA': {'>test_seq1': [(1, 3)]},'CAG': {'>test_seq1': [(2, 4)]},
                            'AGG': {'>test_seq1': [(3, 5)]},'GGT': {'>test_seq1': [(4, 6)]},'GTA': {'>test_seq1': [(5, 7)], '>test_seq2': [(0, 2)]},
                            'TAC': {'>test_seq1': [(6, 8)]},'TAA': {'>test_seq2': [(1, 3)]},'AAC': {'>test_seq2': [(2, 4)]}}
       
       self.assertEqual(id_kmer_dict, exp_id_kmer_dict)
       self.assertEqual(id_read_dict, exp_id_read_dict)
       self.assertEqual(id_presuffix_dict, exp_id_presuffix_dict)
       self.assertEqual(nodes_unique, exp_nodes_unique)
       self.assertEqual(edges_all, exp_edges_all)
       self.assertEqual(suffixes, exp_suffixes)
       self.assertEqual(kmer_id_coord, exp_kmer_id_coord)

if __name__ == '__main__':
    unittest.main()       


