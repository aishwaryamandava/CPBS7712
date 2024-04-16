#!/usr/bin/env python
import unittest
import sys
sys.path.append('../src')
from contigs import concat_contig

class TestConcat(unittest.TestCase):

    def test_concat_contig(self):

        visited = ['TTG','TGC','GCT','CTA','TAG','AGT','GTC','TCA','CAT','ATA']

        exp_contig='TTGCTAGTCATA'

        contig=concat_contig(visited)

        self.assertEqual(contig,exp_contig)
    
if __name__ == '__main__':
    unittest.main()