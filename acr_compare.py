import numpy as np
import pandas as pd
try:
    from statistics import NormalDist
except ModuleNotFoundError:
    from scipy.stats import norm as NormalDist
from math import log
import sys

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
    bins = get_full_union(seg1, seg2)

    # call calc_overlap on each "bin", for both alleles
    # can I assume allele marker (major vs. minor) is same for both inputs? todo
    overlap_scores = np.asarray([b.get_overlap() for b in bins]).flatten()
    overlap_weights = np.asarray([b.get_weights() for b in bins]).flatten()
    weighted_average = np.average(overlap_scores, weights=overlap_weights)

    lengths = np.asarray([b.get_lengths() for b in bins]).sum(axis=0)
    length_names = ['overlap', 'unique_1', 'unique_2']
    length_dict = dict(zip(length_names, lengths))

    # return weighted average
    return weighted_average, length_dict


def calc_overlap(mu1, sigma1, mu2, sigma2):
    """
    Calculate the overlap between two normal distributions, defined by given statistics.

    :param mu1: mean of distribution 1
    :param sigma1: std dev of distribution 1
    :param mu2: mean of distribution 2
    :param sigma2: std dev of distribution 2
    :return: overlap between two distributions as value -> [0, 1]
    """
    # run builtin method in statistics.NormalDist
    if 'statistics' in sys.modules:
        return NormalDist(mu1, sigma1).overlap(NormalDist(mu2, sigma2))

    # if distributions are equal
    if mu1 == mu2 and sigma1 == sigma2:
        return 1

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

def get_full_union(seg1_df, seg2_df):
    full_bins = get_union(seg1_df, seg2_df)
    bins_2 = get_union(seg2_df, seg1_df)

    for bin in bins_2:
        if not bin.unique:
            continue
        else:
            bin.switch_stats()
            full_bins.append(bin)
    return full_bins


def get_union(seg1_df: pd.DataFrame, seg2_df: pd.DataFrame):
    """
    Gets the union of the segment partitions between the two files.

    Throws out any boundaries that do not occur within a segment of the other seg file,
    i.e. creates a bin for every segment held in common

    :param seg1_df: Dataframe of the first seg file
    :param seg2_df: Dataframe of the second seg file
    :return: union of two seg files as a list of Bin objects
    """
    bins = []
    for chrom in np.arange(1, 23):  # assumes segments in each chromosome (except for XY)
        segments1 = seg1_df.loc[seg1_df['Chromosome'] == chrom].reset_index(drop=True)
        segments2 = seg2_df.loc[seg2_df['Chromosome'] == chrom].reset_index(drop=True)

        pointer2 = 0  # starting at beginning of seg2

        for i in range(len(segments1)):
            start = segments1.loc[i]['Start.bp']
            end = segments1.loc[i]['End.bp']

            # create bins for this segment1, updating pointer2 to save computation
            bin_list, pointer2 = create_bins(start, end, segments2, pointer2, segments1.loc[i][STAT_COLUMNS], chrom)
            bins += bin_list

    return bins


def create_bins(start, end, segments2, pointer2, segments1_stats, chrom, bin_list=None):
    """
    Creates bins over this segment defined by start to end, as compared to segments2 df.

    :param chrom:
    :param start: starting base pair (of seg1)
    :param end: ending base pair (of seg1)
    :param segments2: dataframe for seg2
    :param pointer2: current index of seg2 to save computation
    :param segments1_stats: statistics of seg1 to copy to Bin
    :param bin_list: list of bins, for recursion
    :return: final list of Bin
    """
    if not bin_list:
        bin_list = []

    # exit statement for end of segment2
    if pointer2 >= len(segments2):
        return bin_list, pointer2

    start2 = segments2.loc[pointer2]['Start.bp']
    end2 = segments2.loc[pointer2]['End.bp']

    if start <= start2:
        if end <= start2:
            this_bin = Bin(start, end, segments1_stats, None, chrom)  # unique segment
            bin_list.append(this_bin)
            return bin_list, pointer2
        else:
            this_bin = Bin(start, start2, segments1_stats, None, chrom)  # unique segment
            bin_list.append(this_bin)
            if end < end2:
                this_bin = Bin(start2, end, segments1_stats, segments2.loc[pointer2][STAT_COLUMNS], chrom)
                bin_list.append(this_bin)
                return bin_list, pointer2
            else:
                this_bin = Bin(start2, end2, segments1_stats, segments2.loc[pointer2][STAT_COLUMNS], chrom)
                bin_list.append(this_bin)
                return create_bins(end2, end, segments2, pointer2 + 1, segments1_stats, chrom, bin_list)
    else:
        if end < end2:
            this_bin = Bin(start, end, segments1_stats, segments2.loc[pointer2][STAT_COLUMNS], chrom)
            bin_list.append(this_bin)
            return bin_list, pointer2
        else:
            if start < end2:
                this_bin = Bin(start, end2, segments1_stats, segments2.loc[pointer2][STAT_COLUMNS], chrom)
                bin_list.append(this_bin)
                return create_bins(end2, end, segments2, pointer2 + 1, segments1_stats, chrom, bin_list)
            else:
                return create_bins(start, end, segments2, pointer2 + 1, segments1_stats, chrom, bin_list)


class Bin:
    """
    Bin object defining the statistics and length for one Bin
    """

    def __init__(self, start, end, stats1, stats2, chromosome):
        self.length = int(end) - int(start)
        self.unique = stats1 is None or stats2 is None
        self.stats1 = stats1
        self.stats2 = stats2
        self.chromosome = chromosome

    def get_overlap(self):
        """
        Calculates the overlap for both alleles based on the statistics

        :return: list of major and minor overlap scores
        """
        if self.unique:
            return [0, 0]
        major = calc_overlap(self.stats1['mu.major'],
                             self.stats1['sigma.major'],
                             self.stats2['mu.major'],
                             self.stats2['sigma.major'])
        minor = calc_overlap(self.stats1['mu.minor'],
                             self.stats1['sigma.minor'],
                             self.stats2['mu.minor'],
                             self.stats2['sigma.minor'])

        return [major, minor]

    def get_weights(self):
        """
        Returns the two weights for this bin

        :return: list of length weights
        """
        return [self.length] * 2

    def get_lengths(self):
        """
        Returns length in appropriate position of list depending on bin type for summing.

        :return: tuple of (overlap_length, unique_1_length, unique_2_length)
        """
        lengths = [0, 0, 0]
        if not self.unique:
            lengths[0] = self.length
        elif self.stats2 is None:
            lengths[1] = self.length
        elif self.stats1 is None:
            lengths[2] = self.length

        return lengths

    def switch_stats(self):
        store_stat = self.stats1
        self.stats1 = self.stats2
        self.stats2 = store_stat
