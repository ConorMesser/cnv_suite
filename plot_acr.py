import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np


def plot_acr(seg1, seg2, bins, seg1_name, seg2_name):
    min_ylim = 1.2

    fig, ax = plt.subplots(3, figsize=(20, 16), gridspec_kw={'height_ratios': [25, 25, 4]})
    seg_dfs = [seg1, seg2]
    combined_chrom_ends = pd.DataFrame()

    for i in range(2):
        Ag = seg_dfs[i].groupby("Chromosome")
        chr_ends = Ag["End.bp"].apply(max)
        chr_ends[0] = 0
        chr_ends = chr_ends.sort_index()

        combined_chrom_ends = combined_chrom_ends.append(chr_ends, ignore_index=True)

    combined_chrom_ends = combined_chrom_ends.max(axis=0).cumsum()

    for i in range(2):
        Ag = seg_dfs[i].groupby("Chromosome")
        for chrom, C in Ag:
            for _, x in C.iterrows():
                span = x[["Start.bp", "End.bp"]] + combined_chrom_ends[chrom - 1]

                # plot mus
                ax[i].plot(
                  span,
                  x["mu.major"]*np.r_[1, 1],
                  color = 'r',
                  linestyle = (0, [1, 1]),
                  linewidth = 3
                )
                ax[i].plot(
                  span,
                  x["mu.minor"]*np.r_[1, 1],
                  color = 'b',
                  linestyle = (1, [1, 1]),
                  linewidth = 3
                )
                # draw solid line across each dotted segment for clarity
                ax[i].plot(
                  span,
                  x["mu.major"]*np.r_[1, 1],
                  color = 'r',
                  linewidth = 0.5
                )
                ax[i].plot(
                  span,
                  x["mu.minor"]*np.r_[1, 1],
                  color = 'b',
                  linewidth = 0.5
                )

                # plot sigmas
                ax[i].add_patch(patches.Rectangle(
                  [span[0], x["mu.major"] - x["sigma.major"]],
                  span[1] - span[0], 2*x["sigma.major"],
                  color = "r",
                  alpha = 0.5,
                  linewidth = 0
                ))
                ax[i].add_patch(patches.Rectangle(
                  [span[0], x["mu.minor"] - x["sigma.minor"]],
                  span[1] - span[0], 2*x["sigma.minor"],
                  color = "b",
                  alpha = 0.5,
                  linewidth = 0
                ))

            ax[i].axvline(combined_chrom_ends[chrom], linestyle = ":", color = 'k')

            ax[i].set_xlim([0, combined_chrom_ends.iloc[-1]])
            # ax[i].set_ylim(ensure_min_ylim(ax[i].get_ylim(), min_ylim))

    ax[0].set_xticks([])
    ax[1].set_xticks((combined_chrom_ends.iloc[1:].values + combined_chrom_ends[:-1].values) / 2)
    ax[1].set_xticklabels(combined_chrom_ends.iloc[1:].index)
    ax[0].tick_params(labelsize=10)
    ax[1].tick_params(labelsize=10)

    ax[0].set_ylabel(seg1_name)
    ax[1].set_ylabel(seg2_name)

    scores = []
    x_heatmap_mat = []
    bin_line_artists = []
    prev_end = -100000000
    # sort bins by chromosome and start.bp
    bins.sort_values(['chromosome', 'Start.bp'], inplace=True)
    bins.reset_index(inplace=True)

    # grey dashed line for each bin (how to notate which file when non-overlapping: color?)
    for ind in bins.index:  # vectorize todo
        bin_chrom = bins['chromosome'][ind]
        bin_start = bins['Start.bp'][ind] + combined_chrom_ends[bin_chrom - 1]
        bin_end = bins['End.bp'][ind] + combined_chrom_ends[bin_chrom - 1]

        if bins['unique'][ind]:
            if bins['length_1_unique'][ind] > 0:
                color = "#3D3D8C"
                zorder = 0.5
            else:
                color = "#C94643"
                zorder = 0.5
        else:
            color = "black"
            zorder = 0.6

        bin_line_artists.append(ax[1].axvline(bin_end, ymin=0, ymax=2.055, c=color, clip_on=False, linewidth=0.5, alpha=0.5, zorder=zorder))

        if (bin_start - prev_end) > 100:
            bin_line_artists.append(ax[1].axvline(bin_start, ymin=0, ymax=2.055, c=color, clip_on=False, linewidth=0.5, alpha=0.5, zorder=zorder))

            x_heatmap_mat.append(bin_start)
            if not prev_end < 0:
                scores.append(-1)

        x_heatmap_mat.append(bin_end)
        score = (bins['major_overlap'][ind] + bins['minor_overlap'][ind]) / 2  # average of two allelic scores
        scores.append(score)

        prev_end = bin_end

    # heatmap with overlap score for each bin
    y_heatmap_mat = [[0] * len(x_heatmap_mat), [1] * len(x_heatmap_mat)]
    x_heatmap_mat = [x_heatmap_mat, x_heatmap_mat]
    scores = [scores]

    ax[2].pcolormesh(x_heatmap_mat, y_heatmap_mat, scores, cmap='RdGy')  # 2-D irregularly sized array
    ax[2].set_xticks([])
    ax[2].set_yticks([])
    # ax[2].set_frame_on = False
    ax[2].set_ylabel('Comparison Score')

    ax[0].set_ylim(ensure_min_ylim(ax[0].get_ylim(), min_ylim))
    ax[1].set_ylim(ensure_min_ylim(ax[1].get_ylim(), min_ylim))
    # sm = plt.cm.ScalarMappable(cmap='binary', norm=plt.Normalize(vmin=0, vmax=1))
    # sm._A = []
    # plt.colorbar(sm, ax=ax[2])
    fig.tight_layout()

    fig.savefig('full_plot.png')

    ninety_quartile_1 = bins['mu.major_1'].quantile(0.9, interpolation='lower')
    ninety_quartile_2 = bins['mu.major_2'].quantile(0.9, interpolation='lower')
    ax[0].set_ylim(ensure_min_ylim([0, ninety_quartile_1], min_ylim))
    ax[1].set_ylim(ensure_min_ylim([0, ninety_quartile_2], min_ylim))

    bin_line_artists = np.asarray(bin_line_artists)
    change_to_gray = lambda x: x.set_color('gray')
    np.vectorize(change_to_gray)(bin_line_artists)

    plt.show()
    fig.savefig('zoom_plot.png')


def ensure_min_ylim(lim, min_lim):
    if lim[0] > 0:
        lim = (0, lim[1])
    if lim[1] < min_lim:
        return lim[0], min_lim
    else:
        return lim
