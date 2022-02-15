# cnv-methods
Various methods for dealing with copy number data, primarily through visualization, comparison, and simulation.

Comparison: Compare the allelic copy ratio profiles from two segment file, outputting an overlap score and plots. Allows for comparison of the copy number pipeline within different Whole-Exome or Whole-Genome pipelines. Primary usage depends on acr_compare function that takes in two file names.

Simulation: Generate a random copy number profile based on phylogenetic history and overlapping CN events. Apply profile to normal samples (SNVs VCF and coverage .bed files) to generate reference sample with known copy number events for testing and comparing pipelines/tools/technologies.

Visualization: Produce either static or interactive plots to display the CN profile (of a single sample or multiple).
