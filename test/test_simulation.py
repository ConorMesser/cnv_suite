import unittest

import pandas as pd
import numpy as np

from cnv_suite.simulate import CNV_Profile


class MyTestCase(unittest.TestCase):

    def setUp(self) -> None:
        csize = {'1': 50000, '2': 100000}
        cent_loc = {'1': 25000, '2': 30000}
        self.cnv_profile = CNV_Profile(num_subclones=1, csize=csize, cent_loc=cent_loc)

        self.cnv_profile_full = CNV_Profile(num_subclones=1, csize=csize, cent_loc=cent_loc)
        # set phylogeny

        self.cnv_profile_full.phylogeny.ccfs = {1: 1, 2: 0.5}

        self.cnv_profile_full.add_arm(1, p_whole=1, p_q=1, chrom='1', p_deletion=0, allele='maternal')
        self.cnv_profile_full.add_focal(1, chrom='1', p_deletion=0, allele='maternal', position=(1, 2000), cnv_level=5)
        self.cnv_profile_full.add_arm(1, p_whole=0, p_q=1, chrom='2', p_deletion=0, allele='paternal')
        self.cnv_profile_full.add_focal(2, chrom='2', p_deletion=1, allele='maternal', position=(1, 5000), cnv_level=1)

    # test a few of the add_ methods (deterministic)
    def test_add_focal(self):

        self.assertEqual(len(self.cnv_profile.event_trees.keys()), 2)
        self.assertEqual(self.cnv_profile.event_trees['1'].maternal_tree[100].pop().data.cn_change, 1)
        self.assertEqual(self.cnv_profile.event_trees['1'].paternal_tree[100].pop().data.cn_change, 1)

        # test amplification
        self.cnv_profile.add_focal(1, chrom='1', p_deletion=0, allele='maternal', position=(100, 1000), cnv_level=3)
        self.assertEqual(len(self.cnv_profile.event_trees.keys()), 2)
        maternal_set, p_iter = get_event(self.cnv_profile.event_trees['1'].maternal_tree[100], 'focal')
        self.assertEqual(maternal_set.data.cn_change, 3)
        self.assertEqual(self.cnv_profile.event_trees['1'].paternal_tree[100].pop().data.cn_change, 1)

        # test deletion
        self.cnv_profile.add_focal(1, chrom='2', p_deletion=1, allele='paternal', position=(1000, 2000), cnv_level=1)
        self.assertEqual(self.cnv_profile.event_trees['2'].maternal_tree[1000].pop().data.cn_change, 1)
        paternal_set, p_iter = get_event(self.cnv_profile.event_trees['2'].paternal_tree[1000], 'focal')
        self.assertEqual(paternal_set.data.cn_change, -1)

    def test_add_arm(self):
        # test amplification
        self.cnv_profile.add_arm(2, p_whole=0, p_q=1, chrom='1', p_deletion=0, allele='maternal')
        maternal_p, mp_iter = get_event(self.cnv_profile.event_trees['1'].maternal_tree[10000], 'arm')
        maternal_q, mq_iter = get_event(self.cnv_profile.event_trees['1'].maternal_tree[30000], 'arm')
        paternal_p, pp_iter = get_event(self.cnv_profile.event_trees['1'].paternal_tree[10000], 'arm')
        paternal_q, pq_iter = get_event(self.cnv_profile.event_trees['1'].paternal_tree[30000], 'arm')

        self.assertEqual(maternal_p.data.cn_change, 1)
        self.assertEqual(mq_iter, 0)
        self.assertEqual(pp_iter, 0)
        self.assertEqual(pq_iter, 0)

        # test deletion
        self.cnv_profile.add_arm(2, p_whole=1, p_q=1, chrom='2', p_deletion=1, allele='paternal')
        maternal_p, mp_iter = get_event(self.cnv_profile.event_trees['2'].maternal_tree[10000], 'arm')
        maternal_q, mq_iter = get_event(self.cnv_profile.event_trees['2'].maternal_tree[40000], 'arm')
        paternal_p, pp_iter = get_event(self.cnv_profile.event_trees['2'].paternal_tree[10000], 'arm')
        paternal_q, pq_iter = get_event(self.cnv_profile.event_trees['2'].paternal_tree[40000], 'arm')

        self.assertEqual(mp_iter, 0)
        self.assertEqual(mq_iter, 0)
        self.assertEqual(pp_iter, 1)
        self.assertEqual(pq_iter, 1)
        self.assertEqual(paternal_p.data.cn_change, -1)
        self.assertEqual(paternal_q.data.cn_change, -1)

    # def test_add_wgd(self):
    #     self.cnv_profile.add_focal(1, chrom='1', p_deletion=0, allele='maternal', position=(100, 1000), cnv_level=3)
    #     self.cnv_profile.add_arm(2, p_whole=1, p_q=0, chrom='2', p_deletion=0, allele='maternal')
    #
    #     # test amplification
    #     self.cnv_profile.add_wgd(1, both_alleles=True)
    #     maternal_1f, mp_iter = get_event(self.cnv_profile.event_trees['1'].maternal_tree[500], 'arm')
    #     maternal_1a, mq_iter = get_event(self.cnv_profile.event_trees['1'].maternal_tree[2000], 'arm')
    #     maternal_2a, pp_iter = get_event(self.cnv_profile.event_trees['2'].maternal_tree[10000], 'arm')
    #
    #     self.assertEqual(maternal_1f.data.cn_change, 4)
    #     self.assertEqual(maternal_1a.data.cn_change, 1)
    #     self.assertEqual(maternal_2a.data.cn_change, 2)

    # test gen_coverage
    def test_calculate_profiles(self):
        self.cnv_profile_full.calculate_profiles()

        # test dfs
        true_df = pd.DataFrame([['1', 1, 2000, 7, 1],
                                ['1', 2000, 50000, 2, 1],
                                ['2', 1, 5000, 2, 0.5],
                                ['2', 5000, 30000, 2, 1],
                                ['2', 30000, 100000, 1, 1]], columns=['Chromosome', 'Start.bp', 'End.bp', 'mu.major', 'mu.minor'])
        np.testing.assert_equal(self.cnv_profile_full.cnv_profile_df.values, true_df.values)


def get_event(interval_set, event_type):
    event = None
    type_iter = 0
    for p in interval_set:
        if p.data.type == event_type:
            type_iter += 1
            event = p
    return event, type_iter


if __name__ == '__main__':
    unittest.main()
