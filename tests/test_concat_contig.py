#!/usr/bin/env python
import unittest
import sys
sys.path.append('../src')
from contigs import concat_contig

class TestConcat(unittest.TestCase):

    def test_concat_contig(self):

        visited=['AT', 'TC', 'CA', 'AG', 'GG', 'GT', 'TA', 'AA', 'AC']
        exp_contig='ATCAGGTAAC'

        contig=concat_contig(visited)

        self.assertEqual(contig,exp_contig)
    
if __name__ == '__main__':
    unittest.main()