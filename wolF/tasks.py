import wolf


class SimulateProfile(wolf.Task):
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


class CallDepth(wolf.Task):
    inputs = {
        "bam": None,
        "bai": None,
        "bed": None,
        "extra_flags": "-s -a"
    }
    script = """
    samtools view -M -L ${bed} -u ${bam} | samtools depth -b ${bed} ${extra_flags} -o coverage.bed - 
    """
    output_patterns = {"coverage": "coverage.bed"}
    docker = "gcr.io/broad-getzlab-workflows/base_image:v0.0.5"


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

    PROPORTION=`awk -v a="$desired_coverage" -v b="$COV" 'BEGIN {printf a/b }' </dev/null`
    
    if awk -v prop="$PROPORTION" 'BEGIN { exit (prop < 1.0) }' /dev/null
    then  # inverted exit conditions (prop >= 1.0)
        echo "Desired Coverage of ${desired_coverage}x is $PROPORTION times the current bam coverage. Using original bam."
        samtools view $bam -h -o downsample.bam
    else  # prop < 1.0
        INT_FRAC=`awk -v prop="$PROPORTION" -v seed="$subsample_seed" 'BEGIN {printf prop + seed }' </dev/null`
        samtools view $bam -h -s $INT_FRAC -o downsample.bam  # -s INT.FRAC where INT is seed, FRAC is proportion of reads
    fi
    
    samtools index downsample.bam downsample.bam.bai
    """
    output_patterns = {"bam_downsample": "downsample.bam",
                       "bai_downsample": "downsample.bam.bai"}
    docker = "gcr.io/broad-getzlab-workflows/base_image:v0.0.5"
