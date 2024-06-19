import unittest
from cnv_suite.visualize import plot_acr_static, plot_acr_interactive, \
    update_cnv_scatter_sigma_toggle, update_cnv_scatter_color, update_cnv_scatter_cn, make_cnv_scatter
import pandas as pd
import numpy as np
from matplotlib.testing.compare import compare_images
import matplotlib.pyplot as plt


class TestVisualize(unittest.TestCase):

    def setUp(self) -> None:
        csize = {'1': 249250621, '2': 243199373, '3': 198022430, '4': 191154276, '5': 180915260,
                 '6': 171115067, '7': 159138663, '8': 146364022, '9': 141213431, '10': 135534747,
                 '11': 135006516, '12': 133851895, '13': 115169878, '14': 107349540, '15': 102531392,
                 '16': 90354753, '17': 81195210, '18': 78077248, '19': 59128983, '20': 63025520,
                 '21': 48129895, '22': 51304566, '23': 156040895, '24': 57227415}
        self.csize = csize
        seg_df = pd.read_csv('./test/visualize/viz_seg.seg', sep='\t')
        self.seg_df = seg_df

        # make_cnv_scatter(series, fig, lw=0.015, sigmas=False)
        fig, seg_df, trace_start, trace_end = plot_acr_interactive(seg_df, csize,
                         segment_colors='red_blue', sigmas=False,
                         purity=None, ploidy=None,
                         min_seg_lw=0.01, y_upper_lim=2)

        self.plotly_fig = fig

    # test static plotting
    def test_static_plot(self):
        fig, ax = plt.subplots(1, 1, figsize=(20, 5))
        plot_acr_static(self.seg_df, ax, self.csize,
                        segment_colors='Difference', sigmas=False, min_seg_lw=4, y_upper_lim=2)
        test_fn = './test/visualize/static_fig_test.png'
        plt.savefig(test_fn)
        compare_images('./test/visualize/static_fig.png', test_fn, tol=1.0)

    # test plot_acr_interactive todo
    def test_interactive_plot(self):
        pass

    # test make_cnv_scatter
    def test_cnv_scatter(self):
        ser = pd.Series(['1', 1, 5000, 3, 2, 0.5, 1, 5000], ['Chromosome', 'genome_start', 'genome_end',
                                           'mu_major', 'mu_minor', 'sigma_major', 'Start.bp', 'End.bp'])
        make_cnv_scatter(ser, self.plotly_fig, lw=0.02, sigmas=False)

        np.testing.assert_almost_equal(self.plotly_fig.data[-4]['y'], (2, 2, 2, 2))
        np.testing.assert_almost_equal(self.plotly_fig.data[-3]['y'], (3, 3, 3, 3))
        np.testing.assert_almost_equal(self.plotly_fig.data[-2]['y'], (2.02, 1.98, 1.98, 2.02))
        np.testing.assert_almost_equal(self.plotly_fig.data[-1]['y'], (3.02, 2.98, 2.98, 3.02))

        ser = pd.Series(['1', 1, 5000, 3, 2, 0.5, 1, 5000], ['Chromosome', 'genome_start', 'genome_end',
                                           'mu_major', 'mu_minor', 'sigma_major', 'Start.bp', 'End.bp'])
        make_cnv_scatter(ser, self.plotly_fig, lw=0.02, sigmas=True)

        np.testing.assert_almost_equal(self.plotly_fig.data[-4]['y'], (2.5, 1.5, 1.5, 2.5))
        np.testing.assert_almost_equal(self.plotly_fig.data[-3]['y'], (3.5, 2.5, 2.5, 3.5))
        np.testing.assert_almost_equal(self.plotly_fig.data[-2]['y'], (2.02, 1.98, 1.98, 2.02))
        np.testing.assert_almost_equal(self.plotly_fig.data[-1]['y'], (3.02, 2.98, 2.98, 3.02))


    # test update methods
    # update_cnv_scatter_sigma_toggle
    def test_sigma_toggle(self):
        sigmas = self.plotly_fig.select_traces(selector={'name': 'cnv_sigma'})
        for s in sigmas:
            self.assertFalse(s['visible'])
        update_cnv_scatter_sigma_toggle(self.plotly_fig, True)

        for s in sigmas:
            self.assertTrue(s['visible'])

    # update_cnv_scatter_color
    def test_update_color(self):
        self.assertEqual(self.plotly_fig.data[12]['fillcolor'], '#2C38A8')
        self.assertEqual(self.plotly_fig.data[13]['fillcolor'], '#E6393F')
        self.assertEqual(self.plotly_fig.data[14]['fillcolor'], '#2C38A8')
        self.assertEqual(self.plotly_fig.data[15]['fillcolor'], '#E6393F')
        self.assertEqual(self.plotly_fig.data[16]['fillcolor'], '#2C38A8')
        self.assertEqual(self.plotly_fig.data[17]['fillcolor'], '#E6393F')

        update_cnv_scatter_color(self.plotly_fig, ['green'], ['yellow'], 12, 16)
        self.assertEqual(self.plotly_fig.data[12]['fillcolor'], 'green')
        self.assertEqual(self.plotly_fig.data[13]['fillcolor'], 'yellow')
        self.assertEqual(self.plotly_fig.data[14]['fillcolor'], 'green')
        self.assertEqual(self.plotly_fig.data[15]['fillcolor'], 'yellow')
        self.assertEqual(self.plotly_fig.data[16]['fillcolor'], '#2C38A8')
        self.assertEqual(self.plotly_fig.data[17]['fillcolor'], '#E6393F')

    # update_cnv_scatter_cn
    def test_update_cn(self):
        self.assertEqual(self.plotly_fig.data[12]['y'], (0.6, 0.4, 0.4, 0.6))
        self.assertEqual(self.plotly_fig.data[13]['y'], (1.6, 1.4, 1.4, 1.6))
        self.assertEqual(self.plotly_fig.data[14]['y'], (0.51, 0.49, 0.49, 0.51))
        self.assertEqual(self.plotly_fig.data[15]['y'], (1.51, 1.49, 1.49, 1.51))

        update_cnv_scatter_cn(self.plotly_fig, [2.5], [0.3], [0.2], 12, 16, lw=0.02)
        np.testing.assert_almost_equal(self.plotly_fig.data[12]['y'], (0.5, 0.1, 0.1, 0.5))
        np.testing.assert_almost_equal(self.plotly_fig.data[13]['y'], (2.7, 2.3, 2.3, 2.7))
        np.testing.assert_almost_equal(self.plotly_fig.data[14]['y'], (0.32, 0.28, 0.28, 0.32))
        np.testing.assert_almost_equal(self.plotly_fig.data[15]['y'], (2.52, 2.48, 2.48, 2.52))


if __name__ == '__main__':
    unittest.main()
