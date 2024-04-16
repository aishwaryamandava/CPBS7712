#!/usr/bin/env python
import unittest
import sys
sys.path.append('../src')
from alignment import edit_distance

class TestEditDistance(unittest.TestCase):

    def test_edit_distance(self):
        
        query_seq = 'CAGGAAC'
        contig_1 = 'ATCAGGTAAC'
        contig_2 = 'GTAAC'

        exp_edist_1=1
        exp_edist_2=1

        edist_1=edit_distance(query_seq, contig_1)
        edist_2=edit_distance(query_seq, contig_2)

        self.assertEqual(edist_1, exp_edist_1)
        self.assertEqual(edist_2, exp_edist_2)

if __name__ == '__main__':
    unittest.main()