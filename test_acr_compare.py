import unittest
import pandas as pd
from statistics import NormalDist  # only run on version >= Python 3.8

from acr_compare import calc_overlap, create_bins


class CalcOverlap(unittest.TestCase):

    def test_same(self):
        self.assertEqual(calc_overlap(10, 0.5, 10, 0.5), 1)

    def test_far(self):
        self.assertEqual(calc_overlap(0, 0.5, 1000, 0.5), 0)
        self.assertEqual(calc_overlap(10, 0.0001, 20, 0.0001), 0)

    def test_exact(self):
        self.assertEqual(calc_overlap(5, 1, 10, 1), NormalDist(10, 1).cdf(7.5) * 2)
        self.assertEqual(calc_overlap(5, 0.25, 5.1, 0.2),
                         NormalDist(5, 0.25).overlap(NormalDist(5.1, 0.2)))

    def test_order(self):
        self.assertEqual(calc_overlap(10, 0.5, 11, 0.4), calc_overlap(11, 0.4, 10, 0.5))


class Bins(unittest.TestCase):

    def setUp(self) -> None:
        self.segments2 = pd.DataFrame([[10, 20, 1, 1, 1, 1],
                                       [30, 40, 2, 2, 2, 2],
                                       [41, 45, 3, 3, 3, 3]],
                                      columns=['Start.bp', 'End.bp', 'mu.minor',
                                               'sigma.minor', 'mu.major', 'sigma.major'])
        self.one_stats = pd.Series([1, 1, 1, 1], ['mu.minor', 'sigma.minor', 'mu.major', 'sigma.major'])
        self.two_stats = pd.Series([2, 2, 2, 2], ['mu.minor',
                                                                'sigma.minor', 'mu.major', 'sigma.major'])
        self.three_stats = pd.Series([3, 3, 3, 3], ['mu.minor',
                                                                'sigma.minor', 'mu.major', 'sigma.major'])

        self.five_stats = pd.Series([5, 5, 5, 5], ['mu.minor',
                                               'sigma.minor', 'mu.major', 'sigma.major'])

    def test_empty(self):
        self.assertEqual(create_bins(1, 10, self.segments2, 0, self.five_stats, chrom), ([], 0))
        self.assertEqual(create_bins(20, 30, self.segments2, 0, self.five_stats, chrom), ([], 1))

    def test_single(self):
        Bin, pointer = create_bins(10, 15, self.segments2, 0, self.five_stats, chrom)
        self.assertEqual(len(Bin), 1)
        self.assertEqual(Bin[0].length, 5)
        self.assertTrue(Bin[0].stats1.equals(self.five_stats))
        self.assertTrue(Bin[0].stats2.equals(self.one_stats))
        self.assertEqual(pointer, 0)

        Bin, pointer = create_bins(8, 25, self.segments2, 0, self.five_stats, chrom)
        self.assertEqual(len(Bin), 1)
        self.assertEqual(Bin[0].length, 10)
        self.assertTrue(Bin[0].stats1.equals(self.five_stats))
        self.assertTrue(Bin[0].stats2.equals(self.one_stats))
        self.assertEqual(pointer, 1)

        Bin, pointer = create_bins(32, 40, self.segments2, 0, self.five_stats, chrom)
        self.assertEqual(len(Bin), 1)
        self.assertEqual(Bin[0].length, 8)
        self.assertTrue(Bin[0].stats1.equals(self.five_stats))
        self.assertTrue(Bin[0].stats2.equals(self.two_stats))
        self.assertEqual(pointer, 2)

    def test_multiple(self):
        Bin, pointer = create_bins(18, 42, self.segments2, 0, self.five_stats, chrom)
        self.assertEqual(len(Bin), 3)
        self.assertEqual(Bin[0].length, 2)
        self.assertTrue(Bin[0].stats1.equals(self.five_stats))
        self.assertTrue(Bin[0].stats2.equals(self.one_stats))
        self.assertEqual(Bin[1].length, 10)
        self.assertTrue(Bin[1].stats1.equals(self.five_stats))
        self.assertTrue(Bin[1].stats2.equals(self.two_stats))
        self.assertEqual(Bin[2].length, 1)
        self.assertTrue(Bin[2].stats1.equals(self.five_stats))
        self.assertTrue(Bin[2].stats2.equals(self.three_stats))
        self.assertEqual(pointer, 2)




