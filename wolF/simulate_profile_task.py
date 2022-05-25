#!/usr/bin/env python3

import wolf
import pandas as pd
import numpy as np


class Simulate_Profile(wolf.Task):
    inputs = {
        "cnv_pickle",
        "purity",
        "coverage_file",
        "vcf_file",
        "read_depths"
    }
    script = """
    ./cnv_profile.py ${cnv_pickle} ${coverage_file} ${vcf_file} ${read_depths} ${purity} -oc simulated_coverage.txt -oh simulated_hets.txt
    """
    output_patterns = {"coverage": "simulated_coverage.txt",
                       "hets": "simulated_hets.txt"}
    docker = "gcr.io/broad-getzlab-workflows/simulate_cnv_profile:v0.0.1"
