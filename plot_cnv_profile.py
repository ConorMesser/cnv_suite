import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches


def plot_acr(seg_df, sample_name, csize, segment_colors='difference'):
    fig, ax = plt.subplots()
    ax.set_title(f'{sample_name} ACR Plot')

    seg_df['Start.bp'] = seg_df['Start.bp'].astype(int)
    seg_df['End.bp'] = seg_df['End.bp'].astype(int)

    seg_df, chr_order, chrom_start = prepare_df(seg_df, csize, df_type='seg')
    add_background(ax, chr_order, csize, height=7)

    from matplotlib import colors
    cmap = colors.LinearSegmentedColormap.from_list("", ["blue","purple","red"])

    if segment_colors == 'black':
        seg_df['color_bottom'] = '#000000'
        seg_df['color_top'] = '#000000'
    elif segment_colors == 'difference':
        seg_df['color_bottom'] = seg_df.apply(lambda x: cmap(int(np.floor(max(0, (0.5 - 0.5 * calc_color(x['major'] - x['minor'])) * 255)))), axis=1)
        seg_df['color_top'] = seg_df.apply(lambda x: cmap(int(np.floor(min(255, (0.5 + 0.5 * calc_color(x['major'] - x['minor'])) * 255)))), axis=1)

    ax.hlines(seg_df['minor'].values, seg_df['genome_start'], seg_df['genome_end'], color=seg_df['color_bottom'], lw=5)
    ax.hlines(seg_df['major'].values, seg_df['genome_start'], seg_df['genome_end'], color=seg_df['color_top'], lw=5)

    ax.set_xticks(np.asarray(list(chrom_start.values())[:-1]) + np.asarray(list(csize.values())) / 2)
    ax.set_xticklabels(chr_order, fontsize=12)
    ax.set_xlim(chrom_start[chr_order[0]], chrom_start['Z'])

    ax.set_yticks([0, 1, 2])
    ax.set_yticklabels(['0', '1', '2'], fontsize=12)
    ax.set_ylim(-0.05, 5)
    plt.setp(ax.spines.values(), color='gray', lw=.5)
    plt.setp([ax.get_xticklines(), ax.get_yticklines()], color='gray', lw=.5)
    plt.xlabel("Chromosome")
    plt.ylabel("Allelic Copy Ratio")

    return fig


def add_background(ax, chr_order, csize, height=10**7):
    base_start = 0
    chrom_ticks = []
    patch_color = 'white'
    for chrom in chr_order:
        p = patches.Rectangle((base_start, -0.2), csize[chrom], height, fill=True, facecolor=patch_color, edgecolor=None,
                              alpha=.1)  # Background
        ax.add_patch(p)
        patch_color = 'gray' if patch_color == 'white' else 'white'
        chrom_ticks.append(base_start + csize[chrom] / 2)
        base_start += csize[chrom]


def calc_color(mu_diff):
    return (7*mu_diff**2) / (7*mu_diff**2 + 10)


def prepare_df(df, csize, df_type='maf'):
    chr_order = np.asarray(list(csize.keys()))
    chrom_start = {chrom: start for (chrom, start) in
                   zip(np.append(chr_order, 'Z'), [0] + np.insert(np.cumsum([csize[a] for a in chr_order]), 0, 0))}

    df['Chromosome'] = df['Chromosome'].astype(str)
    df = df[df['Chromosome'].isin(chr_order)]

    if df_type == 'maf':
        suffix = '_position'
    else:
        suffix = '.bp'

    df['genome_start'] = df.apply(lambda x: chrom_start[str(x['Chromosome'])] + x[f'Start{suffix}'], axis=1)
    df['genome_end'] = df.apply(lambda x: chrom_start[str(x['Chromosome'])] + x[f'End{suffix}'], axis=1)

    return df, chr_order, chrom_start

