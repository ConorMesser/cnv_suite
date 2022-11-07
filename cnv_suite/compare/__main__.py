#!/bin/bash/env python3

import pandas as pd
import matplotlib.pyplot as plt

from . import num_segments, breakpoint_distance, \
    compare_length_distribution, mu_sigma_difference, acr_compare, plot_acr_comparison


def main():
    # parse args
    import argparse
    import os

    parser = argparse.ArgumentParser(description='Compare two copy number profiles')
    parser.add_argument("profile_one_fn",
                        help='Filename for first CN Segment Profile (by default, the control profile)')
    parser.add_argument("profile_two_fn", help='Filename for second CN Segment Profile')
    parser.add_argument("output_dir", help='Path for directory to save results')

    parser.add_argument("--num_segments", default=False, action="store_true",
                        help='Compare number of segments in the profiles')
    parser.add_argument("--compare_length_dist", default=False, action="store_true",
                        help='Compare the segment length distributions')
    parser.add_argument("--mu_sigma_diff", default=False, action="store_true",
                        help='Compare the differences in the segment mus (and sigmas)')
    parser.add_argument("--breakpoint_dist", default=False, action="store_true",
                        help='Compare the distances between breakpoints')
    parser.add_argument("--overlap_score", default=False, action="store_true",
                        help='Compare the overlap scores between the profiles')
    parser.add_argument("--all", default=False, action="store_true",
                        help='Run all comparison tools')

    # sample names - optional for compare_length_distribution
    parser.add_argument("-sn", "--sample_names", nargs=2)

    # mu_lim, sigma_lim - optional for mu_sigma_difference
    parser.add_argument("--mu_lim", help='Limit for the y (mu) axis')
    parser.add_argument("--sigma_lim", help='Limit for the x (sigma) axis')

    args = parser.parse_args()

    sample_names = ['Profile_1', 'Profile_2'] if not args.sample_names else args.sample_names

    if not os.path.isdir(args.output_dir):
        os.mkdir(args.output_dir)

    if args.num_segments or args.all:
        num_seg_one = num_segments(file_name=args.profile_one_fn)
        num_seg_two = num_segments(file_name=args.profile_two_fn)

        # save to file
        with open(os.path.join(args.output_dir, 'segment_numbers.txt'), 'w') as f:
            f.write(f'{sample_names[0]}\t{sample_names[1]}\n')
            f.write(f'{num_seg_one}\t{num_seg_two}\n')

    if args.compare_length_dist or args.all:
        length_dist_pvalue, length_dist_fig = compare_length_distribution(file_1=args.profile_one_fn,
                                                                          file_2=args.profile_two_fn,
                                                                          sample_names=args.sample_names)
        # save fig
        plt.savefig(os.path.join(args.output_dir, 'lengths_distribution_fig.svg'))

        # save pvalue
        with open(os.path.join(args.output_dir, 'length_distribution_pval.txt'), 'w') as f:
            f.write(f'{length_dist_pvalue}\n')

    if args.mu_sigma_diff or args.all:
        mu_sigma_fig, mu_sigma_ax = mu_sigma_difference(file_1=args.profile_one_fn,
                                                        file_2=args.profile_two_fn,
                                                        mu_lim=args.mu_lim, sigma_lim=args.sigma_lim)

        # save fig
        plt.savefig(os.path.join(args.output_dir, 'mu_sigma_fig.svg'))

    if args.breakpoint_dist or args.all:
        breakpoint_fig, breakpoint_df = breakpoint_distance(file_control=args.profile_one_fn,
                                                            file_case=args.profile_two_fn)

        # save fig (plotly)
        breakpoint_fig.write_image(os.path.join(args.output_dir, 'breakpoint_distances_fig.svg'))

        # save df
        breakpoint_df.to_csv(os.path.join(args.output_dir, 'breakpoint_distances_counts.txt'), sep='\t')

    if args.overlap_score or args.all:
        overlap_score, optimal_scale_factor, optimal_purity, \
        non_overlap_length, overlap_length, bins = acr_compare(file_1=args.profile_one_fn,
                                                               file_2=args.profile_two_fn)

        fig = plot_acr_comparison(pd.read_csv(args.profile_one_fn, sep='\t'),
                                  pd.read_csv(args.profile_two_fn, sep='\t'),
                                  bins, sample_names[0], sample_names[1], args.output_dir)

        # save to file
        with open(os.path.join(args.output_dir, 'overlap_metrics.txt'), 'w') as f:
            f.write(f'overlap_score\toverlap_length\tnon-overlap_length\n')
            f.write(f'{overlap_score}\t{overlap_length}\t{non_overlap_length}\n')


if __name__ == '__main__':
    main()
