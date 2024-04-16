#!/usr/bin/env python
import unittest
import sys
sys.path.append('../src')
from alignment import kmer_contig

class TestKmerContig(unittest.TestCase):

    def test_kmer_contig(self):

        k_len = 3

        longest_contig = "ATCAGGTAAC"
        kmer_id_coord = {
            "ATC": {">test_seq1": [(0,2)]},
            "TCA": {">test_seq1": [(1,3)]},
            "CAG": {">test_seq1": [(2,4)]},
            "AGG": {">test_seq1": [(3,5)]},
            "GGT": {">test_seq1": [(4,6)]},
            "GTA": {">test_seq1": [(5,7)],">test_seq2":[(0,3)]},
            "TAA": {">test_seq2": [(1,3)]},
            "AAC": {">test_seq2": [(2,4)]}}

        exp_kmer_contig_dict = {
            "ATC": {">test_seq1": {
                    "seq_coords": [(0,2)],
                    "contig_coords": [(0,2)]
                }},
            "TCA": {">test_seq1": {
                    "seq_coords": [(1,3)],
                    "contig_coords": [(1,3)]
                }},
            "CAG": {">test_seq1": {
                    "seq_coords": [(2,4)],
                    "contig_coords": [(2,4)]
                }}, 
            "AGG": {">test_seq1": {
                    "seq_coords": [(3,5)],
                    "contig_coords": [(3,5)]
                }},
            "GGT": {">test_seq1": {
                    "seq_coords": [(4,6)],
                    "contig_coords": [(4,6)]
                }},
            "GTA": {">test_seq1": {
                    "seq_coords": [(5,7)],
                    "contig_coords": [(5,7)]
                },">test_seq2": {
                    "seq_coords": [(0,3)],
                    "contig_coords": [(5,7)]}},
            "TAA": {">test_seq2": {
                    "seq_coords": [(1,3)],
                    "contig_coords": [(6,8)]
                }},
            "AAC": {">test_seq2": {
                    "seq_coords": [(2,4)],
                    "contig_coords": [(7,9)]
                }},
        }

        kmer_contig_dict = kmer_contig(k_len, longest_contig, kmer_id_coord)

        self.assertEqual(kmer_contig_dict, exp_kmer_contig_dict)

if __name__ == '__main__':
    unittest.main()