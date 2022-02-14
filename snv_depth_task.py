#!/usr/bin/env python3

import wolf
import pandas as pd
import numpy as np


class CallDepth(wolf.Task):
    inputs = {
        "bam": None,
        "bai": None,
        "bed": None,
        "index_start": 1,
        "index_end": 40000000000,  # 40 billion is greater than number of bases in whole genome
        "extra_flags": "-s -a"
    }
    script = """
    true_end=$(( `wc -l < ${bed} | xargs` < ${index_end} ? `wc -l < ${bed} | xargs` : ${index_end} ))
    sed -n "${index_start},${true_end}p" ${bed} > scatter.bed
    samtools depth -b scatter.bed ${extra_flags} ${bam} > coverage.bed
    """
    output_patterns = {"coverage": "coverage.bed"}
    docker = "gcr.io/broad-getzlab-workflows/base_image:v0.0.5"


def depth_workflow(bam, bai, bed, scatter_num=100, with_localization=True):
    if with_localization:
        localization = wolf.localization.LocalizeToDisk(
            files={
                "bam": bam,
                "bai": bai,
                "bed": bed
            }
        )
    else:
        localization = dict("bam": bam, "bai": bai, "bed", bed)
    
    # split bed file into scatter_num even partitions
    bed_df = pd.read_csv(bed, sep='\t')
    line_num = bed_df.shape[0]
    start_indices = np.linspace(1, line_num+1, num=min(line_num, scatter_num)+1, dtype=int)
    end_indices = start_indices - 1
    start_indices = start_indices[:-1].tolist()
    end_indices = end_indices[1:].tolist()
    
    # could also call split_intervals task but probably unnecessary because each interval is only 1 bp long

    snv_depth = CallDepth(
        name="get_read_depths",
        inputs={
            "bam": localization['bam'],
            "bai": localization['bai'],
            "bed": localization['bed'],
            "index_start": start_indices,
            "index_end": end_indices
        }
    )

def coverage_workflow(bam, bai, scatter_num=100, bed_path=None, interval_size=None, coverage_resources=None):
    assert bed_path or interval_size
    
    tumor_bam_localization_task = wolf.localization.LocalizeToDisk(
      files = dict(
        bam = bam,
        bai = bai
      )
    )

    tumor_bam_localization = tumor_bam_localization_task.run()
    
    # split_intervals not working on WES target bed file...need to partition ad hoc
    if bed_path:
        bed_df = pd.read_csv(bed_path, sep='\t', header=None, names=['chr', 'start', 'end'], dtype={0: str, 1: int, 2: int})
        bed_df['len'] = bed_df['end'] - bed_df['start']
        scatter_num = 200
        import numpy as np
        num_partitions_chr  = np.ceil(bed_df.groupby('chr')['len'].sum() / bed_df['len'].sum() * scatter_num)
        intervals_out = []
        for c, i in num_partitions_chr.items():
            this_chr_bed = bed_df.loc[bed_df['chr'] == c].reset_index()
            approx_partition_len = this_chr_bed['len'].sum() / i

            start = this_chr_bed.iloc[0]['start']
            len_cum = 0
            for num, row in this_chr_bed.iterrows():
                len_cum += row['len']
                if num == this_chr_bed.shape[0] - 1:
                    intervals_out.append([c, start, row['end'] + 1])
                elif len_cum > approx_partition_len:
                    intervals_out.append([c, start, row['end'] + 1])
                    start = this_chr_bed.iloc[num+1]['start']
                    len_cum = 0

        subset_intervals = pd.DataFrame(intervals_out, columns=['chr', 'start', 'end'])

    elif interval_size:
        # split target list
        split_intervals = wolf.ImportTask(
          task_path = "git@github.com:getzlab/split_intervals_TOOL.git",
          task_name = "split_intervals"
        )

        split_intervals_task = split_intervals.split_intervals(
          bam = tumor_bam_localization["bam"],
          bai = tumor_bam_localization["bai"],
          interval_type = "bed",
          N = 100,
          selected_chrs = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']
        )

        split_intervals_results = split_intervals_task.run()

        # shim task to transform split_intervals files into subset parameters for covcollect task

        def interval_gather(interval_files):
            ints = []
            for f in interval_files:
                ints.append(pd.read_csv(f, sep = "\t", header = None, names = ["chr", "start", "end"]))
            return pd.concat(ints).sort_values(["chr", "start", "end"])

        subset_intervals = interval_gather(split_intervals_results["interval_files"])

    # get coverage
    cov_collect = wolf.ImportTask(
      task_path = "git@github.com:getzlab/covcollect.git",
      task_name = "covcollect"
    )

    cov_collect_task = cov_collect.Covcollect(
      inputs = dict(
        bam = tumor_bam_localization["bam"],
        bai = tumor_bam_localization["bai"],
        intervals = bed_path if bed_path else interval_size,
        subset_chr = subset_intervals["chr"],
        subset_start = subset_intervals["start"],
        subset_end = subset_intervals["end"],
      ),
      retry=2,
      resources=coverage_resources
    )

    cov_collect_results = cov_collect_task.run()

    # gather coverage
    cov_gather = wolf.Task(
      name = "gather_coverage",
      inputs = { "coverage_beds" : [cov_collect_results["coverage"]] },
      script = """cat $(cat ${coverage_beds}) > coverage_cat.bed""",
      outputs = { "coverage" : "coverage_cat.bed" }
    )

    cov_gather_results = cov_gather.run()
    
