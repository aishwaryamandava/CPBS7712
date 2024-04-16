#!/usr/bin/env python
import unittest
import sys
sys.path.append('../src')
from contigs import depthfirst_eulerian_path_

class TestDFS(unittest.TestCase):

    def test_depthfirst_eulerian_path_(self):

        graph= {'AT': ['TC'], 'TC': ['CA'],
              'CA': ['AG'],'AG': ['GG'],
              'GG': ['GT'],'GT': ['TA', 'TA'],
              'TA': ['AC', 'AA'],'AA': ['AC']}

        start_1='AT'
        start_2='GT'

        k_len=3

        exp_paths_1 = ['ATCAGGTAAC']
        exp_paths_2 = ['GTAAC']

        paths_1 = depthfirst_eulerian_path_(graph, start_1)
        paths_2 = depthfirst_eulerian_path_(graph, start_2)

        self.assertEqual(paths_1, exp_paths_1)
        self.assertEqual(paths_2, exp_paths_2)

if __name__ == '__main__':
    unittest.main()