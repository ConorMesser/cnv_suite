import pandas as pd
import scipy.stats as ss
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import math

from cnv_helper_methods import get_segment_interval_trees


def num_segments(file_name=None, seg_df=None):
    assert file_name or seg_df is not None

    if file_name:
        return sum(1 for i in open(file_name, 'rb')) - 1
    else:
        return seg_df.shape[0]


def compare_length_distribution(file_1=None, file_2=None, seg_df_1=None, seg_df_2=None):
    assert (file_1 or seg_df_1 is not None) and (file_2 or seg_df_2 is not None)

    if file_1:
        seg_df_1 = pd.read_csv(file_1, sep='\t')
    if file_2:
        seg_df_2 = pd.read_csv(file_2, sep='\t')

    seg_df_1['length'] = seg_df_1['End.bp'] - seg_df_1['Start.bp']
    seg_df_2['length'] = seg_df_2['End.bp'] - seg_df_2['Start.bp']

    mannwhitney = ss.mannwhitneyu(seg_df_1['length'], seg_df_2['length'])

    fig, ax = plt.subplots(1, 2)
    round_to = 10**7
    bins = np.linspace(0, round_to * math.ceil((np.concatenate([seg_df_1['length'].values, seg_df_2['length'].values])).max() / round_to), 50)
    ax[0].hist(seg_df_1['length'], bins, alpha=0.5, label='x', density=True)
    ax[0].hist(seg_df_2['length'], bins, alpha=0.5, label='y', density=True)
    ax[0].legend(['First profile', 'Second profile'], loc='upper right')
    ax[0].set_xlabel('Segment Length (bp)')
    ax[0].set_ylabel('Density of Counts')
    ax[0].set_title('Segment Lengths', fontsize=10)

    # add plot that shows segment length, normalized to chromosome (or arm?)
    chr_size = pd.concat([seg_df_1.groupby('Chromosome')['End.bp'].max(),
                          seg_df_2.groupby('Chromosome')['End.bp'].max()], axis=1)
    chr_size = chr_size.max(axis=1).to_dict()
    seg_df_1['length_normalized'] = seg_df_1.apply(lambda x: x['length'] / chr_size[x['Chromosome']], axis=1)
    seg_df_2['length_normalized'] = seg_df_2.apply(lambda x: x['length'] / chr_size[x['Chromosome']], axis=1)
    bins = np.linspace(0, 1, 50)
    ax[1].hist(seg_df_1['length_normalized'], bins, alpha=0.5, label='x', density=True)
    ax[1].hist(seg_df_2['length_normalized'], bins, alpha=0.5, label='y', density=True)
    ax[1].legend(['First profile', 'Second profile'], loc='upper right')
    ax[1].set_xlabel('Normalized Segment Length')
    ax[1].set_title('Lengths Normalized by Contig Length', fontsize=10)

    fig.suptitle('Comparing Segment Length Distributions', fontsize=16)

    return mannwhitney.pvalue, fig


def breakpoint_distance(file_control=None, file_case=None, seg_df_control=None, seg_df_case=None):
    assert (file_control or seg_df_control is not None) and (file_case or seg_df_case is not None)

    if file_control:
        seg_df_control = pd.read_csv(file_control, sep='\t')
    if file_case:
        seg_df_case = pd.read_csv(file_case, sep='\t')

    # breakpoints can span large distances (over centromere but also within chromosome)
    # should I calculate distance for both sides of breakpoint?
    # does direction matter?

    # visualizations:
    # - can make scatter plot of distance (y axis) over genome
    # - could make scatter plot of distance by some sort of CN difference representation


def mu_sigma_difference(file_1=None, file_2=None, seg_df_1=None, seg_df_2=None, mu_lim=None, sigma_lim=None):
    assert (file_1 or seg_df_1 is not None) and (file_2 or seg_df_2 is not None)

    if file_1:
        seg_df_1 = pd.read_csv(file_1, sep='\t')
    if file_2:
        seg_df_2 = pd.read_csv(file_2, sep='\t')

    if 'Sample_ID' not in seg_df_1:
        seg_df_1['Sample_ID'] = 'profile_1'
    if 'Sample_ID' not in seg_df_2:
        seg_df_2['Sample_ID'] = 'profile_2'

    # should I take the intersection of the segments??
    # I think so, because the larger segment intersections will be prioritized (in terms of their lengths)
    contig_trees = get_segment_interval_trees(pd.concat([seg_df_1, seg_df_2]))

    # get differences, with non-overlapping segments removed (this can be visualized better using acr_compare tool)
    diff_df = []
    for tree in contig_trees:
        diff_df.append(get_differences_from_intervals(tree, sample_names=[seg_df_1.loc[0, 'Sample_ID'],
                                                                          seg_df_2.loc[0, 'Sample_ID']],
                                                      remove_unmatched=True))

    diff_df = pd.concat(diff_df)

    fig, ax = plt.subplots(1,1)

    # code taken from Claudia's cgaprojects_ibm_tAML_analysis
    ax.axvline(0, c='k', linewidth=0.5)
    ax.axhline(0, c='k', linewidth=0.5)

    sigma_diffs = np.concatenate([diff_df['sigma_minor_diff'], diff_df['sigma_major_diff']])
    mu_diffs = np.concatenate([diff_df['mu_minor_diff'], diff_df['mu_major_diff']])
    lengths = np.concatenate([diff_df['length'], diff_df['length']])
    pcm = ax.scatter(sigma_diffs, mu_diffs,
                     c=['blue']*diff_df.shape[0] + ['red']*diff_df.shape[0], alpha=0.2)

    max_length = np.max(lengths)
    pcm.set_sizes(lengths / max_length * 100 + 3)

    max_mu_diff = np.amax(np.abs(mu_diffs))
    max_sigma_diff = np.amax(np.abs(sigma_diffs))
    if mu_lim is not None and sigma_lim is not None:
        ax.set_ylim(-mu_lim, mu_lim)
        ax.set_xlim(-sigma_lim, sigma_lim)
    else:
        ax.set_ylim(-max_mu_diff, max_mu_diff)
        ax.set_xlim(-max_sigma_diff if np.amin(sigma_diffs) < 0 else 0, max_sigma_diff)
    ax.set_ylabel('ACR (mu) Difference')
    ax.set_xlabel('Variance (sigma) Difference')
    ax.set_title('Comparison of CNV Profiles')

    # legend
    size_list = np.arange(0, max_length, np.ceil(max_length / 5))[1:]
    handles = [
        plt.scatter([-1], [-1], s=length / max_length * 100 + 3, c='gray')
        for length
        in size_list]
    labels = ["{:.2e}bp".format(length) for length in size_list]
    ax.legend(handles, labels, loc=3, framealpha=1, frameon=True, title='Segment Length')
    # fig.show()

    return fig, ax


def get_differences_from_intervals(contig_tree, sample_names=None, remove_unmatched=True):
    if sample_names is None:
        sample_names = ['profile_1', 'profile_2']
    length_arr = []
    mu_minor_diff = []
    mu_major_diff = []
    sigma_minor_diff = []
    sigma_major_diff = []

    for interval in contig_tree:
        if len(interval.data) == 1:
            if remove_unmatched:
                continue
            else:
                d = list(interval.data.values())[0]  # we know there is only one value

                if sample_names[0] in interval.data.keys():
                    mu_minor_diff.append(1 - d.mu_minor)
                    mu_major_diff.append(1 - d.mu_major)
                    sigma_minor_diff.append(-d.sigma_minor)
                    sigma_major_diff.append(-d.sigma_major)
                else:
                    mu_minor_diff.append(d.mu_minor - 1)
                    mu_major_diff.append(d.mu_major - 1)
                    sigma_minor_diff.append(d.sigma_minor)
                    sigma_major_diff.append(d.sigma_major)
        elif len(interval.data) == 2:
            mu_minor_diff.append(interval.data[sample_names[1]].mu_minor -
                                 interval.data[sample_names[0]].mu_minor)
            mu_major_diff.append(interval.data[sample_names[1]].mu_major -
                                 interval.data[sample_names[0]].mu_major)
            sigma_minor_diff.append(interval.data[sample_names[1]].sigma_minor -
                                    interval.data[sample_names[0]].sigma_minor)
            sigma_major_diff.append(interval.data[sample_names[1]].sigma_major -
                                    interval.data[sample_names[0]].sigma_major)
        else:
            raise ValueError(f'Interval should have one or two samples, not {len(interval.data)}')

        length_arr.append(interval.end - interval.begin)

    return pd.DataFrame(np.asarray([length_arr, mu_minor_diff, mu_major_diff, sigma_minor_diff, sigma_major_diff]).transpose(),
                        columns=['length', 'mu_minor_diff', 'mu_major_diff', 'sigma_minor_diff', 'sigma_major_diff'])



