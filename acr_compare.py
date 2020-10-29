import numpy as np
import pandas as pd
try:
    from statistics import NormalDist
except ModuleNotFoundError:
    from scipy.stats import norm as NormalDist
from math import log
import sys

from plot_acr import plot_acr

STAT_COLUMNS = ['mu.minor', 'sigma.minor', 'mu.major', 'sigma.major']


def acr_compare(file_1=None, file_2=None):
    """
    Compare the allelic copy ratio segments between the two seg files.

    :param file_1: name and path of file_1 (default to user input)
    :param file_2: name and path of file_2 (default to user input)
    :return: weighted average of all compared ACR segment normal distributions
    """
    # read in two seg files
    if not file_1 or not file_2:
        file_1 = input('Input first file name: ')
        file_2 = input('Input second file name: ')
    seg1 = pd.read_csv(file_1, sep='\t')
    seg2 = pd.read_csv(file_2, sep='\t')

    # format dataframes (throwing out rows with nan values for mu/sigma)
    seg1 = seg1.dropna(subset=STAT_COLUMNS).reset_index(drop=True)
    seg2 = seg2.dropna(subset=STAT_COLUMNS).reset_index(drop=True)

    # take union of segments to allocate bins
    bins = get_union(seg1, seg2)

    # call calc_overlap on each "bin", for both alleles
    # can I assume allele marker (major vs. minor) is same for both inputs? todo
    bins["major_overlap"] = bins.apply(lambda x: calc_overlap(x["mu.major_1"], x["sigma.major_1"],
                                                              x["mu.major_2"], x["sigma.major_2"],
                                                              x['unique']), axis=1)
    bins["minor_overlap"] = bins.apply(lambda x: calc_overlap(x["mu.minor_1"], x["sigma.minor_1"],
                                                              x["mu.minor_2"], x["sigma.minor_2"],
                                                              x['unique']), axis=1)

    # get average overlap score
    overlap_score = np.average(np.concatenate([bins['major_overlap'], bins['minor_overlap']]),
                               weights=np.concatenate([bins['length'], bins['length']]))

    non_overlap_length = int(bins['length_1_unique'].sum()) + int(bins['length_2_unique'].sum())
    overlap_length = int(bins['length_overlap'].sum())

    # plot files with bins marked
    plot_acr(seg1, seg2, bins, 'REBC Seg File', 'wolF Seg File')

    # overlap_scores = np.asarray([b.get_overlap() for b in bins]).flatten()
    # overlap_weights = np.asarray([b.get_weights() for b in bins]).flatten()
    # weighted_average = np.average(overlap_scores, weights=overlap_weights)
    #
    # lengths = np.asarray([b.get_lengths() for b in bins]).sum(axis=0)
    # length_names = ['overlap', 'unique_1', 'unique_2']
    # length_dict = dict(zip(length_names, lengths))

    # return weighted average
    return overlap_score, non_overlap_length, overlap_length, bins

def calc_overlap(mu1, sigma1, mu2, sigma2, unique):
    """
    Calculate the overlap between two normal distributions, defined by given statistics.

    :param mu1: mean of distribution 1
    :param sigma1: std dev of distribution 1
    :param mu2: mean of distribution 2
    :param sigma2: std dev of distribution 2
    :return: overlap between two distributions as value -> [0, 1]
    """
    # check if any values are None (non-overlapping)
    if unique:
        return 0

    # check if distributions are equal
    if mu1 == mu2 and sigma1 == sigma2:
        return 1

    # run builtin method in statistics.NormalDist
    if 'statistics' in sys.modules:
        return NormalDist(mu1, sigma1).overlap(NormalDist(mu2, sigma2))

    # calculate intersection(s) of the two distributions
    x_intersect = calc_pdf_intersect(mu1, sigma1, mu2, sigma2)

    if len(x_intersect) == 1:  # sigma1 == sigma2
        area = NormalDist(mu1, sigma1).cdf(x_intersect[0])
        if area < 0.5:
            return area * 2  # doubled -> pdf1_area = pdf2_area when sigma1 = sigma2
        else:  # take other side of cdf
            return (1 - area) * 2

    # calculate overlap cdf
    mid_section1 = NormalDist(mu1, sigma1).cdf(max(x_intersect)) - \
                   NormalDist(mu1, sigma1).cdf(min(x_intersect))
    mid_section2 = NormalDist(mu2, sigma2).cdf(max(x_intersect)) - \
                   NormalDist(mu2, sigma2).cdf(min(x_intersect))

    # compute sum of overlap sections based on which middle section is larger
    if mid_section1 < mid_section2:
        sum_overlap = 1 + mid_section1 - mid_section2
    else:
        sum_overlap = 1 + mid_section2 - mid_section1

    return sum_overlap


def calc_pdf_intersect(mu1, sigma1, mu2, sigma2):
    """
    Calculate intersection(s) of two normal distributions.

    Returns one value (in a list) if sigmas are equal; otherwise returns two values.

    :param mu1: mean of distribution 1
    :param sigma1: std dev of distribution 1
    :param mu2: mean of distribution 2
    :param sigma2: std dev of distribution 2
    :return: list of intersection values
    """
    if sigma1 == sigma2:
        return [(mu1 + mu2) / 2]

    # calculate roots of log[pdf_1] = log[pdf_2]
    a = 0.5 * (1/sigma2**2 - 1/sigma1**2)
    b = mu1 / sigma1**2 - mu2 / sigma2**2
    c = 0.5 * (mu2**2 / sigma2**2 - mu1**2 / sigma1**2) + log(sigma2 / sigma1)

    return np.roots([a, b, c])


####################

def get_union(seg1_df, seg2_df):  # todo either initialize df with large size or create from list
    full_bins = _union_one_sided(seg1_df, seg2_df)
    bins_2 = _union_one_sided(seg2_df, seg1_df)

    # take only unique values from bins_2
    unique_2_df = bins_2.loc[bins_2['unique'] == True]

    # swap stats columns (and 1/2 length columns)
    l = []
    for c in unique_2_df.columns.values:
        if '_1' in c:
            c = c.replace('1', '2')
        else:
            c = c.replace('2', '1')
        l.append(c)
    unique_2_df.columns = l

    bins = pd.concat([full_bins, unique_2_df], ignore_index=True)

    return bins


def _union_one_sided(seg1_df: pd.DataFrame, seg2_df: pd.DataFrame):
    """
    Gets the union of the segment partitions between the two files.

    Throws out any boundaries that do not occur within a segment of the other seg file,
    i.e. creates a bin for every segment held in common

    :param seg1_df: Dataframe of the first seg file
    :param seg2_df: Dataframe of the second seg file
    :return: union of two seg files as a list of Bin objects
    """
    bin_list = []
    for chrom in np.arange(1, 23):  # assumes segments in each chromosome (except for XY)
        segments1 = seg1_df.loc[seg1_df['Chromosome'] == chrom].reset_index(drop=True)
        segments2 = seg2_df.loc[seg2_df['Chromosome'] == chrom].reset_index(drop=True)

        pointer2 = 0  # starting at beginning of seg2

        for i in range(len(segments1)):
            start = segments1.loc[i]['Start.bp']
            end = segments1.loc[i]['End.bp']

            # create bins for this segment1, updating pointer2 to save computation
            bin_list, pointer2 = create_bins(start, end, segments2, pointer2, segments1.loc[i][STAT_COLUMNS],
                                             chrom, bin_list=bin_list)

    bin_df = pd.DataFrame(bin_list)
    return bin_df


def create_bins(start, end, segments2, pointer2, segments1_stats, chrom, bin_list=None):
    """
    Creates bins over this segment defined by start to end, as compared to segments2 df.

    :param chrom: this chromosome
    :param start: starting base pair (of seg1)
    :param end: ending base pair (of seg1)
    :param segments2: dataframe for seg2
    :param pointer2: current index of seg2 to save computation
    :param segments1_stats: statistics of seg1 to copy to Bin
    :param bin_list: list of bins, for recursion
    :return: final list of Bin
    """
    if bin_list is None:
        bin_list = []

    # exit statement for end of segment2
    if pointer2 >= len(segments2):
        return bin_list, pointer2

    start2 = segments2.loc[pointer2]['Start.bp']
    end2 = segments2.loc[pointer2]['End.bp']

    if start <= start2:
        if end <= start2:
            bin_list.append(append_bin(start, end, segments1_stats, None, chrom))  # unique segment
            return bin_list, pointer2
        else:
            bin_list.append(append_bin(start, start2, segments1_stats, None, chrom))  # unique segment

            if end < end2:
                bin_list.append(append_bin(start2, end, segments1_stats, segments2.loc[pointer2][STAT_COLUMNS], chrom))
                return bin_list, pointer2
            else:
                bin_list.append(append_bin(start2, end2, segments1_stats, segments2.loc[pointer2][STAT_COLUMNS], chrom))
                return create_bins(end2, end, segments2, pointer2 + 1, segments1_stats, chrom, bin_list=bin_list)
    else:
        if end < end2:
            bin_list.append(append_bin(start, end, segments1_stats, segments2.loc[pointer2][STAT_COLUMNS], chrom))
            return bin_list, pointer2
        else:
            if start < end2:
                bin_list.append(append_bin(start, end2, segments1_stats, segments2.loc[pointer2][STAT_COLUMNS], chrom))
                return create_bins(end2, end, segments2, pointer2 + 1, segments1_stats, chrom, bin_list=bin_list)
            else:
                return create_bins(start, end, segments2, pointer2 + 1, segments1_stats, chrom, bin_list=bin_list)


def append_bin(start, end, stats1, stats2, chromosome):
    unique = stats1 is None or stats2 is None
    length = end - start
    length_overlap = 0
    length_1_unique = 0
    length_2_unique = 0

    if not unique:
        length_overlap = length
    elif stats2 is None:
        length_1_unique = length
    elif stats1 is None:
        length_2_unique = length

    bin_dict = {'chromosome': chromosome,
                'Start.bp': start,
                'End.bp': end,
                'unique': unique,
                'length': length,
                'length_overlap': length_overlap,
                'length_1_unique': length_1_unique,
                'length_2_unique': length_2_unique}

    for key in STAT_COLUMNS:
        if stats1 is not None:
            val = stats1[key]
        else:
            val = None
        bin_dict[f'{key}_1'] = val

        if stats2 is not None:
            val = stats2[key]
        else:
            val = None
        bin_dict[f'{key}_2'] = val

    # major_overlap, minor_overlap = get_overlap(stats1, stats2, unique)
    # bin_dict['major_overlap'] = major_overlap
    # bin_dict['minor_overlap'] = minor_overlap
    #
    # modDf = df.append(bin_dict, ignore_index=True)  #todo super inefficient

    return bin_dict


def get_overlap(stats1, stats2, unique):
    """
    Calculates the overlap for both alleles based on the statistics

    :return: list of major and minor overlap scores
    """
    if unique:
        return 0, 0
    major = calc_overlap(stats1['mu.major'],
                         stats1['sigma.major'],
                         stats2['mu.major'],
                         stats2['sigma.major'])
    minor = calc_overlap(stats1['mu.minor'],
                         stats1['sigma.minor'],
                         stats2['mu.minor'],
                         stats2['sigma.minor'])

    return major, minor
