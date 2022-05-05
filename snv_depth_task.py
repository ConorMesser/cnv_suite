#!/usr/bin/env python3

import wolf
import pandas as pd
import numpy as np
import subprocess
import prefect

from wolf.localization import LocalizeToDisk, DeleteDisk


# simulate_cnv_profile = wolf.ImportTask(
#   task_path = "git@github.com:getzlab/cnv_methods_TOOL.git",
#   task_name = "simulate_cnv_profile"
# )


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


def depth_workflow(bed, scatter_num=100, bam=None, bai=None, tumor_bam_localization=None):
    assert (bam and bai) or tumor_bam_localization
    
    if not tumor_bam_localization:
        tumor_bam_localization = wolf.localization.LocalizeToDisk(
            files={
                "bam": bam,
                "bai": bai
            })
    
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
            "bam": tumor_bam_localization['bam'],
            "bai": tumor_bam_localization['bai'],
            "bed": bed,
            "index_start": start_indices,
            "index_end": end_indices
        }
    )
    
    ### todo filter/take subset based on exome bait intervals
    
    if not tumor_bam_localization:
        DeleteDisk(
          inputs = {
            "disk" : [tumor_bam_localization["bam"]],
            "upstream" : snv_depth["coverage"]
          }
        )
    
    return snv_depth


def coverage_workflow(scatter_num=100, bam=None, bai=None, bed_path=None, interval_size=None, coverage_resources=None, tumor_bam_localization=None, genome_file=None):
    assert bed_path or interval_size
    assert (bam and bai) or tumor_bam_localization

    if not tumor_bam_localization:
        tumor_bam_localization = wolf.localization.LocalizeToDisk(
          files = dict(
            bam = bam,
            bai = bai
          ))
    
    # Partition the padded WES intervals list into even parts
    if bed_path and genome_file:
        
        #@prefect.task
        def partition_wes_intervals(bed_path, genome_file, scatter_num):
            bed_df = pd.read_csv(bed_path, sep='\t', header=None, names=['chr', 'start', 'end'], dtype={0: str, 1: int, 2: int})
            bed_df['len'] = bed_df['end'] - bed_df['start']
            scatter_num = scatter_num
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
                        
            return pd.DataFrame(intervals_out, columns=['chr', 'start', 'end'])

        subset_intervals = partition_wes_intervals(bed_path, genome_file, scatter_num)

    else:  # for wgs or wes if genome file not provided
        split_intervals = wolf.ImportTask(
          task_path = "git@github.com:getzlab/split_intervals_TOOL.git",
          task_name = "split_intervals"
        )

        split_intervals_task = split_intervals.split_intervals(
          bam = tumor_bam_localization["bam"],
          bai = tumor_bam_localization["bai"],
          interval_type = "bed",
          N = scatter_num,
          selected_chrs = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']
        )

        # shim task to transform split_intervals files into subset parameters for covcollect task
        #@prefect.task
        def interval_gather(interval_files):
            ints = []
            for f in interval_files:
                ints.append(pd.read_csv(f, sep = "\t", header = None, names = ["chr", "start", "end"]))
            return pd.concat(ints).sort_values(["chr", "start", "end"])

        subset_intervals = interval_gather(split_intervals_task["interval_files"])

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

    # gather coverage
    cov_gather = wolf.Task(
      name = "gather_coverage",
      inputs = { "coverage_beds" : [cov_collect_task["coverage"]] },
      script = """cat $(cat ${coverage_beds}) > coverage_cat.bed""",
      outputs = { "coverage" : "coverage_cat.bed" }
    )
    
    if not tumor_bam_localization:
        DeleteDisk(
            inputs = {
            "disk" : [tumor_bam_localization["bam"]],
            "upstream" : cov_collect_task["coverage"]
            }
        )
    
    return cov_gather
    
    
class DownSampleBam(wolf.Task):
    inputs = {
        "bam": None,
        "desired_coverage": None,
        "exome_bed": "",
        "coverage_size": "",
        "subsample_seed": "1"
    }
    script = """
    # get read length
    READ_LEN=$(samtools view ${bam} | head -n100 | cut -f6 | \
    perl -ne '@ciglens = split(/[^0-9]+/); $tot = 0;
    foreach(@ciglens) { $tot += $_; }
    $maxtot = $tot > $maxtot ? $tot : $maxtot;
    END { print $maxtot }')
    
    if [[ -n $exome_bed ]]
    then
        COV_SIZE=$(awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' ${exome_bed})
    elif [[ -n $coverage_size ]]
    then
        COV_SIZE=$coverage_size
    else
        COV_SIZE=""
    fi
    
    if [[ -z $COV_SIZE ]]
    then
        # get approximate mean coverage - over whole genome
        COV=$(samtools idxstats ${bam} | awk -v read_len=$READ_LEN '
        NR <= 24 { genome_len += $2; n_reads += $3 }
        END { print read_len*n_reads/genome_len }')
    else
        # get approx mean coverage - over specified coverage area
        COV=$(samtools idxstats ${bam} | awk -v read_len=$READ_LEN -v cov_size=$COV_SIZE '
        NR <= 24 { n_reads += $3 }
        END { print read_len*n_reads/cov_size }')
    fi

    PROPORTION=$(echo "scale=4 ; $desired_coverage / $COV" | bc)
    
    if [ `echo "$PROPORTION > 1.0" | bc` == 1 ]
    then
        echo "Desired Coverage of ${desired_coverage}x is $PROPORTION times the current bam coverage. Using original bam."
        cat $bam > downsample.bam
    else
        samtools view $bam -h -s (echo $prop + $seed | bc) > downsample.bam  # -s INT.FRAC where INT is seed, FRAC is proportion of reads
    fi
    
    samtools index downsample.bam downsample.bam.bai
    """
    output_patterns = {"bam_downsample": "downsample.bam",
                       "bai_downsample": "downsample.bam.bai"}
    docker = "gcr.io/broad-getzlab-workflows/base_image:v0.0.5"
    
    
def full_simulation_workflow(bam, bai, vcf_bed, cnv_pickle, purity,
                             scatter_num_cov=100, scatter_num_depth=100, 
                             wes_target_intervals=None, interval_size=None, 
                             coverage_resources=None, desired_coverage=None):
    #
    # Downsample bam to specified level of coverage
    if desired_coverage:
        downsample_bam = DownSampleBam(name="downsample_bam",
            inputs={
                "bam": bam,
                "desired_coverage": desired_coverage,
                "exome_bed": wes_target_intervals if wes_target_intervals else ""
            })
    
    #
    # Localize bam and bai files, using either given bam or downsampled bam
    tumor_bam_localization = LocalizeToDisk(
      files = dict(
        bam = downsample_bam["bam_downsample"] if desired_coverage else bam,
        bai = downsample_bam["bai_downsample"] if desired_coverage else bai
      ))
    
    #
    # Run Coverage workflow, calling cov_collect on either WGS intervals or a WES target bed file
    coverage_results = coverage_workflow(scatter_num=scatter_num_cov, bed_path=wes_target_intervals, interval_size=interval_size, coverage_resources=coverage_resources, tumor_bam_localization=tumor_bam_localization)
    
    # Shim task to make bed file from input VCF
    @prefect.task
    def vcf_to_bed(input_vcf):
        full_vcf_df = pd.read_csv('../HapASeg/simulation/NA12878/VCF/NA12878.vcf', sep='\t', comment='#', header=None, 
                         names=['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','NAME'])
        bed_df = full_vcf_df[['CHROM', 'POS']]
        bed_df['end_pos'] = bed_df['POS']
        bed_df['POS'] = bed_df['POS'] - 1  # bed is zero-based, but samtools in one-based (output depth will include end but not start)
        output_fn = './vcf.bed'
        bed_df.to_csv(output_fn, sep='\t', index=False, header=False)
        return output_fn
    
    #
    # Run SNV Depth workflow, calling samtools depth on every site in given VCF bed file
    snv_depth_results = depth_workflow(vcf_to_bed(input_vcf), scatter_num=scatter_num_depth, tumor_bam_localization=tumor_bam_localization)        
    
    
    #
    # Create CNV Profile based on given cnv_pickle and save Coverage and VCF files for given purity
    cnv_profile_results = simulate_cnv_profile.Simulate_Profile(
        inputs = dict(
            cnv_pickle = cnv_pickle,
            purity = purity,
            coverage_file = coverage_results["coverage"],
            vcf_file = input_vcf,
            read_depths = snv_depth_results["coverage"]
        )
    )
    
    #
    # Run HapASeg pipeline
    #
    hapaseg_results = hapaseg_workflow.workflow(
  callstats_file = None,

  tumor_coverage_bed = None,
  normal_coverage_bed = None,

  ref_genome_build='hg38', #must be hg19 or hg38
  
  target_list = None,
  localization_token='',

  num_cov_seg_samples=5,

  persistant_dry_run = False,
        VCF_file = input_vcf,
        tumor_allele_counts = cnv_profile_results["coverage"],
        # todo get normal
        # normal_allele_counts = cnv_profile_results["coverage"]
        
        # need to add in hapaseg_load_snps_task["allele_counts"] rather than using callstats/mutect
)
    
    #
    # Compare results to true copy number profile
    
    
    
    DeleteDisk(
        inputs = {
        "disk" : [tumor_bam_localization["bam"]],
        "upstream" : cov_collect_task["coverage"]
        }
    )
    
    
    
    
    
    
    
 
