#!/usr/bin/env python
import unittest
import sys
sys.path.append('../src')
from alignment import seed_extend_align

class TestSeedExtendAlign(unittest.TestCase):

    def test_seed_extend_align(self):

        query_seq = 'CAGGAAC'
        seed_len = 3

        dfs_contig_list = {'contig_1': 'ATCAGGTAAC',
                           'contig_2': 'GTAAC'}
        
        query_to_kmer_index={'CAG': [0], 'GAA': [3]}

        exp_align_query = {
            'contig1':{'index':[[0]],'edist':1,'contig_len':10,'edist/contig_len':0.1},
            'contig2':{'index':[],'edist':0,'contig_len':5,'edist/contig_len':0.2}
        }

        exp_sorted_align_query={}

        exp_longest_contig=None
        
        sorted_align_query, longest_contig = seed_extend_align(query_seq, seed_len, dfs_contig_list, query_to_kmer_index)

        self.assertEqual(sorted_align_query, exp_sorted_align_query)
        self.assertEqual(longest_contig, exp_longest_contig)

if __name__ == '__main__':
    unittest.main()

