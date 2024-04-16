#!/usr/bin/env python
import unittest
import sys
sys.path.append('../src')
from alignment import aln_file

class TestAlnFile(unittest.TestCase):

    def test_aln_file(self):

        kmer_contig_dict = {
            "ATC": {">test_seq1": {"seq_coords": [(0,2)], "contig_coords": [(0,2)]}},
            "TCA": {">test_seq1": {"seq_coords": [(1,3)], "contig_coords": [(1,3)]}},
            "CAG": {">test_seq1": {"seq_coords": [(2,4)], "contig_coords": [(2,4)]}}, 
            "AGG": {">test_seq1": {"seq_coords": [(3,5)], "contig_coords": [(3,5)]}},
            "GGT": {">test_seq1": {"seq_coords": [(4,6)], "contig_coords": [(4,6)]}},
            "GTA": {
                ">test_seq1": {"seq_coords": [(5,7)], "contig_coords": [(5,7)]},
                ">test_seq2": {"seq_coords": [(0,2)], "contig_coords": [(5,7)]}
            },
            "TAA": {">test_seq2": {"seq_coords": [(1,3)], "contig_coords": [(6,8)]}},
            "AAC": {">test_seq2": {"seq_coords": [(2,4)], "contig_coords": [(7,9)]}}
        }
        sorted_align_query = {'contig_1': {'index': [[0]],'edist': 1,
                                           'contig_len': 10,'edist/contig_len': 0.1},
                            'contig_2': {'index': [], 'edist': 0, 
                                         'contig_len': 5, 'edist/contig_len': 0}}


        exp_aln_list = [('test_seq2', 'contig_1', 0, 4, 5, 9), ('test_seq1', 'contig_1', 0, 7, 0, 7)]

        aln_list = aln_file(kmer_contig_dict, sorted_align_query)

        self.assertEqual(aln_list, exp_aln_list)

if __name__=='__main__':
    unittest.main()
