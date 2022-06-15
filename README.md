# CNV_Suite
Various methods for dealing with copy number data, primarily through visualization, comparison, and simulation.

- Comparison: Compare the allelic copy ratio profiles from two segment files. Allows for comparison of the copy number pipeline within different Whole-Exome or Whole-Genome pipelines or using different CN tools.

- Simulation: Generate a random copy number profile based on phylogenetic history and overlapping CN events. Apply profile to normal samples (SNVs VCF and coverage .bed files) to generate reference sample with known copy number events for testing and comparing pipelines/tools/technologies.

- Visualization: Produce static or interactive plots to display the CN profile (of a single sample or multiple). 

## Installation

Clone this repo into your environment: `git clone git@github.com:getzlab/cnv_suite.git`

Install from the setup.py script: `pip install .`

Requirements include *pandas*, *numpy*, *scipy*, *matplotlib*, *plotly*, *intervaltree*, *natsort*, and *pandarallel*.

## Command Line Usage

### Compare
Runs various comparison tools on the two given ACR segment profiles. Choose `--all` to run all tools.
```
compare First_profile_filename Second_profile_filename [--num_segments] [--compare_length_dist] [--mu_sigma_diff] [--breakpoint_dist] [--all] [--sample_names] [--mu_lim] [--sigma_lim]
```
For more details, run `compare -h`

### Simulate
Simulate the VCF read depths and binned coverage according to the given criteria. Requires a pickle file of a CNV_profile object created using the `cnv_suite.simulate` package.
```
simulate cnv_pickle_file coverage_file vcf_file read_depths purity [--output_coverage filename] [--output_hets filename] [--normal_coverage coverage_file] [--normal_depths read_depths]
```
For more details, run `simulate -h`

### Visualize
Save a static CNV plot for the given segment profile with the given options.
```
visualize segment_filename output_filename [--csize_file] [--segment_colors] [--hide_sigmas] [--min_seg_lw] [--y_upper]
```
For more details, run `visualize -h`

## Package Usage

Beyond the static plotting, the visualization package allows for the production of interactive plots using plotly (`cnv_suite.visualize.plot_acr_interactive`). Methods also exist for updating the segment colors, copy number values, and sigma visibility (for use in a dashboard or notebook). 

The simulation package provides the `CNV_Profile` class which can generate random (or specified) copy number events (arm level, focal level, chromothripsis, WGD) with haplotype and phylogeny awareness. This profile can then be applied to reference samples to return simulated SNV read counts and binned coverage. 

Finally, other utilities exist to calculate absolute copy number values given purity and ploidy (`cnv_suite.calc_absolute_cn`), the average weighted copy number for a profile (`cnv_suite.calc_avg_cn`), copy number values at a particular loci (`cnv_suite.return_seg_data_at_loci` and `cnv_suite.apply_segment_data_to_df` to apply this to a maf dataframe), and many others. 
