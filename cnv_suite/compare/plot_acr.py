import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import os


def plot_acr_comparison(seg1, seg2, bins, sample_one_name, sample_two_name, output_dir):
    """
    Make ACR comparison plot with both ACR segment plots, bin lines and heatmap for overlap scores.

    :param seg1: first seg dataframe
    :param seg2: second seg dataframe
    :param bins: bin dataframe
    :param sample_one_name: name for first sample
    :param sample_two_name: name for second sample
    :param output_dir: path to file output directory
    :return: ACR comparison figure
    """
    # constant parameters
    min_ylim = 1.2

    # initialize plot
    fig, axes = plt.subplots(3, 2, figsize=(20, 16),
                             gridspec_kw={'height_ratios': [25, 25, 4],
                                          'width_ratios': [75, 1]})
    ax = np.asarray(axes)[:, 0]
    ax_colorbar = np.asarray(axes)[:, 1]
    seg_dfs = [seg1, seg2]

    # sort bins by chromosome and start.bp
    bins.sort_values(['chromosome', 'Start.bp'], inplace=True)
    bins.reset_index(inplace=True)

    # get the greater of each chromosome end
    combined_chr_ends = pd.DataFrame()
    for i in range(2):
        Ag = seg_dfs[i].groupby("Chromosome")
        chr_ends = Ag["End.bp"].apply(max)
        chr_ends[0] = 0
        chr_ends = chr_ends.sort_index()
        combined_chr_ends = combined_chr_ends.append(chr_ends, ignore_index=True)
    combined_chr_ends = combined_chr_ends.max(axis=0).cumsum()

    # plot primary ACR plots
    seg_names = [sample_one_name, sample_two_name]
    for i in range(2):
        if i == 0:
            xticks=[]
        else:
            xticks = (combined_chr_ends.iloc[1:].values + combined_chr_ends[:-1].values) / 2
        plot_acr(seg_dfs[i], ax[i], combined_chr_ends,
                 seg_names[i], xticks)
    # only set xticklabels for bottom plot
    ax[1].set_xticklabels(combined_chr_ends.iloc[1:].index)

    # add vertical line for each bin start/stop and collect scores for heatmap
    scores = []
    x_heatmap_mat = []
    bin_line_artists = []
    prev_end = -1000000
    for chr, start, end, unique, length_1, major, minor in zip(bins['chromosome'], bins['Start.bp'],
                                                               bins['End.bp'], bins['unique'],
                                                               bins['length_1_unique'],
                                                               bins['major_overlap'],
                                                               bins['minor_overlap']):
        bin_start = start + combined_chr_ends[chr - 1]
        bin_end = end + combined_chr_ends[chr - 1]

        if unique:
            zorder = 0.5
            if length_1 > 0:
                color = "#3D3D8C"
            else:
                color = "#C94643"
        else:
            color = "black"
            zorder = 0.6

        # draw bin_end line with color based on bin "type"
        bin_line_artists.append(ax[1].axvline(bin_end, ymin=0, ymax=2.055,
                                              c=color, linewidth=0.5, alpha=0.5,
                                              clip_on=False, zorder=zorder))
        # if gap exists between previous bin and this bin:
        #   draw bin_start line and add score of -1 for empty bin
        if (bin_start - prev_end) > 100:
            bin_line_artists.append(ax[1].axvline(bin_start, ymin=0, ymax=2.055,
                                                  c=color, linewidth=0.5, alpha=0.5,
                                                  clip_on=False, zorder=zorder))
            x_heatmap_mat.append(bin_start)
            if not prev_end < 0:  # no empty bin at start needed todo
                scores.append(-1)

        # this bin score and ending point
        x_heatmap_mat.append(bin_end)
        score = (major + minor) / 2  # average of two allelic scores
        scores.append(score)

        prev_end = bin_end

    # heatmap with overlap score for each bin
    ax[2].pcolormesh([x_heatmap_mat, x_heatmap_mat],
                     [[0] * len(x_heatmap_mat), [1] * len(x_heatmap_mat)],
                     [scores], cmap='RdGy')
    # heatmap layout
    ax[2].set_xticks([])
    ax[2].set_yticks([])
    ax[2].set_ylabel('Comparison Score')

    # add colorbar for heatmap (0 to 1) and separate legend for missing bins
    sm = plt.cm.ScalarMappable(cmap='binary', norm=plt.Normalize(vmin=0, vmax=1))
    sm._A = []
    plt.colorbar(sm, cax=ax_colorbar[2])
    missing_patch = patches.Patch(color='#67001f', label='No segment data')
    ax_colorbar[2].legend(handles=[missing_patch], loc='lower left',
                          bbox_to_anchor=(-0.5, -0.4), frameon=False)
    fig.delaxes(ax_colorbar[0])
    fig.delaxes(ax_colorbar[1])
    fig.tight_layout()
    _ensure_min_ylim(ax[0], ax[0].get_ylim(), min_ylim)
    _ensure_min_ylim(ax[1], ax[1].get_ylim(), min_ylim)

    # save full plot with colored bins
    fig.savefig(os.path.join(output_dir, 'acr_comparison_full_color.svg'))

    # save zoomed plot with colored bins
    ninety_quartile_1 = bins['mu.major_1'].quantile(0.9, interpolation='lower')
    ninety_quartile_2 = bins['mu.major_2'].quantile(0.9, interpolation='lower')
    _ensure_min_ylim(ax[0], [0, ninety_quartile_1], min_ylim)
    _ensure_min_ylim(ax[1], [0, ninety_quartile_2], min_ylim)
    fig.savefig(os.path.join(output_dir, 'acr_comparison_zoom_color.svg'))

    # save zoomed plot with gray bins for clarity
    bin_line_artists = np.asarray(bin_line_artists)
    change_to_gray = lambda x: x.set_color('gray')
    np.vectorize(change_to_gray)(bin_line_artists)
    fig.savefig(os.path.join(output_dir, 'acr_comparison_zoom_grey.svg'))

    return fig


def plot_acr(seg_df, ax, combined_chr_ends, seg_name, xticks):
    """
    Plot primary ACR plot with mus and sigma boxes.
    """
    Ag = seg_df.groupby("Chromosome")
    for chrom, C in Ag:
        for _, x in C.iterrows():
            span = x[["Start.bp", "End.bp"]] + combined_chr_ends[chrom - 1]

            # plot mus
            ax.plot(
                span,
                x["mu.major"] * np.r_[1, 1],
                color='r',
                linestyle=(0, [1, 1]),
                linewidth=3
            )
            ax.plot(
                span,
                x["mu.minor"] * np.r_[1, 1],
                color='b',
                linestyle=(1, [1, 1]),
                linewidth=3
            )
            # draw solid line across each dotted segment for clarity
            ax.plot(
                span,
                x["mu.major"] * np.r_[1, 1],
                color='r',
                linewidth=0.5
            )
            ax.plot(
                span,
                x["mu.minor"] * np.r_[1, 1],
                color='b',
                linewidth=0.5
            )

            # plot sigmas
            ax.add_patch(patches.Rectangle(
                [span[0], x["mu.major"] - x["sigma.major"]],
                span[1] - span[0], 2 * x["sigma.major"],
                color="r",
                alpha=0.7,
                linewidth=0
            ))
            ax.add_patch(patches.Rectangle(
                [span[0], x["mu.minor"] - x["sigma.minor"]],
                span[1] - span[0], 2 * x["sigma.minor"],
                color="b",
                alpha=0.7,
                linewidth=0
            ))

        # set layout
        ax.axvline(combined_chr_ends[chrom], linestyle=":", color='k')  # delineate chromosomes
        ax.set_xlim([0, combined_chr_ends.iloc[-1]])
        ax.tick_params(labelsize=10)
        ax.set_ylabel(seg_name)
        ax.set_xticks(xticks)


def _ensure_min_ylim(ax, lim, min_lim):
    """
    Calculates y limits, setting as (min(0, lim[0]), max(min_lim, lim[1]).
    """
    ax.set_ylim(min(0, lim[0]), max(min_lim, lim[1]))
