#!/usr/bin/env python
import unittest
import sys
sys.path.append('../src')
from contigs import depthfirst_eulerian_path_

class TestDFS(unittest.TestCase):

    def test_depthfirst_eulerian_path_(self):

        graph={
            'AT': ['TA'],'TA': ['AG','AT'],
            'AG': ['GT'],'GT': ['TC'],
            'TC': ['CA'],'CA': ['AT'],
            'TT': ['TG'],'TG': ['GC'],
            'GC': ['CT'],'CT': ['TA']
        }

        start='TT'

        k_len=3

        exp_paths = ['TTGCTAG']

        paths = depthfirst_eulerian_path_(graph, start)

        self.assertEqual(paths, exp_paths)

if __name__ == '__main__':
    unittest.main()