import pandas as pd

import wolf
from wolf.localization import LocalizeToDisk, DeleteDisk
import prefect


cnv_sim = wolf.ImportTask(
    task_path=".",
    task_name="cnv_sim"
)


def coverage_workflow(scatter_num=100, bam=None, bai=None, bed_path=None, interval_size=None, coverage_resources=None, tumor_bam_localization=None, genome_file=None):
    assert bed_path or interval_size
    assert (bam and bai) or tumor_bam_localization

    if not tumor_bam_localization:
        tumor_bam_localization = wolf.localization.LocalizeToDisk(
            files = dict(
                bam = bam,
                bai = bai
            ))

    # Partition the genome into
    split_intervals_task = wolf.ImportTask(
        task_path = "git@github.com:getzlab/split_intervals_TOOL.git",
        task_name = "split_intervals"
    )

    split_intervals = split_intervals_task.split_intervals(
        bam = tumor_bam_localization["bam"],
        bai = tumor_bam_localization["bai"],
        interval_type = "bed",
        N = scatter_num,
        selected_chrs = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']
    )

    # shim task to transform split_intervals files into subset parameters for covcollect task
    @prefect.task
    def interval_gather(interval_files):
        ints = []
        for f in interval_files:
            ints.append(pd.read_csv(f, sep = "\t", header = None, names = ["chr", "start", "end"]))
        return pd.concat(ints).sort_values(["chr", "start", "end"])

    subset_intervals = interval_gather(split_intervals["interval_files"])

    # get coverage
    cov_collect_task = wolf.ImportTask(
        task_path = "git@github.com:getzlab/covcollect.git",
        task_name = "covcollect"
    )

    cov_collect = cov_collect_task.Covcollect(
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
        inputs = { "coverage_beds" : [cov_collect["coverage"]] },
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


def full_simulation_workflow(bam, bai, input_vcf, cnv_pickle, purity,
                             scatter_num_cov=100,
                             wes_target_intervals=None, interval_size=None,
                             coverage_resources=None, desired_coverage=None):
    #
    # Downsample bam to specified level of coverage
    if desired_coverage:
        downsample_bam = cnv_sim.DownSampleBam(name="downsample_bam",
                                       inputs={
                                           "bam": bam,
                                           "desired_coverage": desired_coverage,
                                           "exome_bed": wes_target_intervals if wes_target_intervals else ""
                                       })

    #
    # Localize bam and bai files, using either given bam or downsampled bam
    tumor_bam_localization = LocalizeToDisk(
        files=dict(
            bam=downsample_bam["bam_downsample"] if desired_coverage else bam,
            bai=downsample_bam["bai_downsample"] if desired_coverage else bai
        ))

    #
    # Run Coverage workflow, calling cov_collect on either WGS intervals or a WES target bed file
    coverage_results = coverage_workflow(scatter_num=scatter_num_cov, bed_path=wes_target_intervals,
                                         interval_size=interval_size, coverage_resources=coverage_resources,
                                         tumor_bam_localization=tumor_bam_localization)

    # Shim task to make bed file from input VCF
    @prefect.task
    def vcf_to_bed(input_vcf):
        full_vcf_df = pd.read_csv('../HapASeg/simulation/NA12878/VCF/NA12878.vcf', sep='\t', comment='#', header=None,
                                  names=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT',
                                         'NAME'])
        bed_df = full_vcf_df[['CHROM', 'POS']]
        bed_df['end_pos'] = bed_df['POS']
        bed_df['POS'] = bed_df[
                            'POS'] - 1  # bed is zero-based, but samtools in one-based (output depth will include end but not start)
        output_fn = './vcf.bed'
        bed_df.to_csv(output_fn, sep='\t', index=False, header=False)
        return output_fn

    #
    # Run SNV Depth workflow, calling samtools depth on every site in given VCF bed file
    snv_depth_results = cnv_sim.CallDepth(
        name="get_read_depths",
        inputs={
            "bam": tumor_bam_localization['bam'],
            "bai": tumor_bam_localization['bai'],
            "bed": vcf_to_bed(input_vcf),
        }
    )


    # Create CNV Profile based on given cnv_pickle and save Coverage and VCF files for given purity
    cnv_profile_results = cnv_sim.SimulateProfile(
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
    #     hapaseg_results = hapaseg_workflow.workflow(
    #   callstats_file = None,

    #   tumor_coverage_bed = None,
    #   normal_coverage_bed = None,

    #   ref_genome_build='hg38', #must be hg19 or hg38

    #   target_list = None,
    #   localization_token='',

    #   num_cov_seg_samples=5,

    #   persistant_dry_run = False,
    #         VCF_file = input_vcf,
    #         tumor_allele_counts = cnv_profile_results["coverage"],
    #         # todo get normal
    #         # normal_allele_counts = cnv_profile_results["normal_coverage"]
    #         # normal_allele_counts = cnv_profile_results["normal_coverage"]

    #         # need to add in hapaseg_load_snps_task["allele_counts"] rather than using callstats/mutect
    # )

    #
    # Compare results to true copy number profile

    DeleteDisk(
        inputs={
            "disk": [tumor_bam_localization["bam"]],
            "upstream": [coverage_results["coverage"], snv_depth_results["coverage"]]
        }
    )
