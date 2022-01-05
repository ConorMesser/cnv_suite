#!/usr/bin/env python3

import wolf
import pandas as pd
import numpy as np


class CallDepth(wolf.Task):
    inputs = {
        "bam": None,
        "bed": None,
        "index_start": 1,
        "index_end": 40000000000,  # 40 billion is greater than number of bases
        "extra_flags": "-s -a"
    }
    script = """
    true_end=$(( `wc -l < ${bed} | xargs` < ${index_end} ? `wc -l < ${bed} | xargs` : ${index_end} ))
    sed -n '${index_start},${true_end}p' ${bed} > scatter.bed
    samtools depth -b scatter.bed ${extra_flags} ${bam} > coverage.bed
    """
    output_patterns = {"coverage": "coverage.bed"}
    docker = "gcr.io/broad-getzlab-workflows/base_image:v0.0.5"


def depth_workflow(bam, bai, bed, scatter_num=100):
    localization = wolf.localization.BatchLocalDisk(
        files={
            "bam": bam,
            "bai": bai,
            "bed": bed
        }
    )

    bed_df = pd.read_csv(bed, sep='\t')
    line_num = bed_df.shape[0]
    start_indices = np.linspace(1, line_num, num=min(line_num, scatter_num)+1, dtype=int)
    end_indices = start_indices - 1
    start_indices = start_indices[:-1]
    end_indices = end_indices[1:]

    snv_depth = CallDepth(
        name="get_read_depths",
        inputs={
            "bam": localization['bam'],
            "bed": localization['bed'],
            "index_start": start_indices,
            "index_end": end_indices
        }
    )


if __name__ == '__main__':
    # generate bed from vcf file
    vcf_df = pd.read_csv('./simulation/NA12878/VCF/NA12878.vcf', sep='\t', comment='#')
    bed_df = vcf_df[['CHROM', 'POS']]
    bed_df['end_pos'] = bed_df['POS']
    bed_df.to_csv('./local_vcf.bed', sep='\t', index=False, header=False)

    with wolf.Workflow(
            workflow=depth_workflow,
            conf={"clust_frac": 0.2},  # limit to 20% of available nodes, to avoid going over project quota
            common_task_opts={"retry": 5}  # retry every task up to 5 times
    ) as w2:
        w2.run(
            bam="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/CEU/NA12878/exome_alignment/NA12878.alt_bwamem_GRCh38DH.20150826.CEU.exome.cram",
            bai="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/CEU/NA12878/exome_alignment/NA12878.alt_bwamem_GRCh38DH.20150826.CEU.exome.cram.crai",
            bed="./local_vcf.bed",
            scatter_num=100,

            RUN_NAME="test_WES_NA12878"
        )

    C_all = []
    for c in w2.results.loc[(slice(None), "get_read_depths"), ("outputs", "coverage")]:
        C_all.append(pd.read_csv(c, sep="\t", names=["chr", "pos", "depth"]))

    C_all = pd.concat(C_all)
