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
    
       exp_id_kmer_dict = {'>test_seq1': ['ATA','TAG','AGT','GTC','TCA','CAT'], '>test_seq2': ['TTG','TGC','GCT','CTA','TAT']}
       exp_id_read_dict = {'>test_seq1': 'ATAGTCAT', '>test_seq1': 'TTGCTAT'}
       exp_id_presuffix_dict = {'>test_seq1': [['AT', 'TA'], ['TA', 'AG'],['AG','GT'],['GT','TC'],['TC','CA'],['CA','AT']], 
                                '>test_seq2': [['TT', 'TG'],['TG','GC'],['GC','CT'],['CT','TA'],['TA','AT']]}
       exp_nodes_unique = {'AT', 'TA', 'AG', 'GT', 'TC','CA','TT','TG','GC','CT'}
       exp_edges_all = [('AT', 'TA'), ('TA', 'AG'), ('AG', 'GT'), ('GT', 'TC'), ('TC', 'CA'), ('TT', 'TG'), ('TG', 'GC'), ('GC', 'CT'), ('CT', 'TA')]
       exp_suffixes = {'TA', 'AG', 'GT','TC','CA','AT','TG','GC','CT'}
       exp_kmer_id_coord = {'ATA': {'>test_seq1': [(0, 2)]}, 'TAG': {'>test_seq1': [(1, 3)]}, 'AGT': {'>test_seq1': [(2, 4)]},
                                   'GTC': {'>test_seq1': [(3, 5)]}, 'TCA': {'>test_seq1': [(4, 6)]}, 'CAT': {'>test_seq1': [(5, 7)]},
                                   'TTG': {'>test_seq2': [(0, 2)]}, 'TGC': {'>test_seq2': [(1, 3)]}, 'GCT': {'>test_seq2': [(2, 4)]},
                                   'CTA': {'>test_seq2': [(3, 5)]}, 'TAT': {'>test_seq2': [(4, 6)]}}
       
       self.assertEqual(id_kmer_dict, exp_id_kmer_dict)
       self.assertEqual(id_read_dict, exp_id_read_dict)
       self.assertEqual(id_presuffix_dict, exp_id_presuffix_dict)
       self.assertEqual(nodes_unique, exp_nodes_unique)
       self.assertEqual(edges_all, exp_edges_all)
       self.assertEqual(suffixes, exp_suffixes)
       self.assertEqual(kmer_id_coord, exp_kmer_id_coord)

if __name__ == '__main__':
    unittest.main()       


