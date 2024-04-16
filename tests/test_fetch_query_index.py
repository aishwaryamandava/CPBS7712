#!/usr/bin/env python
import unittest
import sys
sys.path.append('../src')
from alignment import fetch_query_index

class TestFetchIndex(unittest.TestCase):

    def test_fetch_query_index(self):
        query_path='./test_QUERY.fasta'
        seed_len=3

        exp_query_seq='CAGGAAC'
        exp_query_to_kmer_index={'CAG': [0], 'GAA': [3]}

        query_seq, query_to_kmer_index = fetch_query_index(query_path, seed_len)

        self.assertEqual(query_seq, exp_query_seq)
        self.assertEqual(query_to_kmer_index, exp_query_to_kmer_index)

if __name__ == '__main__':
    unittest.main()