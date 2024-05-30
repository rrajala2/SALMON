#!/usr/bin/bash -l
# Disclaimer: This code is for instructional purposes only.
#             You may use and alter this code, but you are responsible for
#               understanding any code you use for research and publications.

# Load the required modules:
module load slurm        # for SLURM scheduler commands (like sbatch)
module load fastqc       # for Quality Control
module load multiqc      # for Quality Control
module load salmon       # for aligning to transcriptome
module load subread      # for featureCounts
module load R            # for DE-analysis
module load bioconductor # for DE-analysis

# Show what modules are loaded (thus capturing versions and previously loaded)
module list

# Bash "strict mode"
set -ue

function usage() {
    cat <<USAGE

    Usage: $0 --fastq-dir <DIR> --ids-file <FILE> --transcripts <FILE> [--is-single (T|F) ] [--run-fastq-qc (T|F) ]

    Options:
        --fastq-dir       Path to the FASTQ dir (e.g. "reads")
        --transcripts     Path to the transcriptome reference FASTA file (e.g. 00-genome/transcripts.fa)
        --ids-file        File containing list of IDs to process
        --minAssignedFrags Minimum fragments required to map to a transcript to be included in the output (defaults to 10)
        --is-single       Indication if this is single-ended reads (otherwise paired-end is assumed)
        --run-fastq-qc    Run FastQC and MultiQC on FASTQ files (defaults to F, use "T" to turn on)
        --reuse-counts    Reuse existing counts files (i.e. don't regenerate counts)
        #TODO: --generate-counts Do everything through count generation steps

    Example:
        $0 --ids-file ids.txt --transcripts 00-genome/transcripts.fa --fastq-dir 00-reads

USAGE
    exit 1
}

threads=12    # Number of CPU cores
mem_per_cpu=6 # 6 GBS of RAM per CPU core
run_fastq_qc="F"
reuse_counts="F"
generate_counts="T"

# Make some extra directories to help file organization:
code_dir=scripts
qc_dir=01-qc
salmon_dir=02-salmon
raw_counts_dir=04-raw_counts
pca_dir=05-pca
de_dir=06-DE_lists
heatmap_dir=07-heatmaps
log_dir=slurm_logs

mkdir -p $qc_dir
mkdir -p $salmon_dir
mkdir -p $raw_counts_dir
mkdir -p $pca_dir
mkdir -p $de_dir
mkdir -p $heatmap_dir
mkdir -p $log_dir


# Step 0 - getting things set up

is_single=''
minAssignedFrags=10 # Currently default in salmon
while [ $# -gt 0 ]; do
    case $1 in              # check for these flags
        --fastq-dir)
            fastq_dir=$2
            ;;
        --ids-file)
            ids_file=$2
            ;;
        --is-single)
            is_single=$( echo $2 | tr '[:lower:]' '[:upper:]' | cut -c1 )
            ;;
        --run-fastq-qc)
            run_fastq_qc=$( echo $2 | tr '[:lower:]' '[:upper:]' | cut -c1 )
            ;;
        --minAssignedFrags)
            minAssignedFrags=$2
            ;;
        --reuse-counts)
            reuse_counts=$( echo $2 | tr '[:lower:]' '[:upper:]' | cut -c1 )
            ;;
        --transcripts)
            transcripts=$2
            ;;
        --help | -h)
            usage
            exit 1
            ;;
        *)
            usage
            exit 1
            ;;
    esac
    # Remove the first flag and value. (i.e. This removes the current $1 and $2,
    #   thus shifting $3 to $1, $4 to $2, etc.)
    shift # Remove the first arg which was the first flag processed
    shift # Remove what is now the first arg, which was the first value processed
done

# Double check that these are set. If not, give the error message.
echo ${transcripts?'Please specify the name of the transcriptome reference file using --transcripts filename (e.g. transcripts.fa)'} > /dev/null
echo ${fastq_dir?'Please specify the FASTQ directory using --fastq-dir dirname'} > /dev/null
echo ${ids_file?'Please specify the name of the IDs file using --ids-file filename (e.g. ids.txt)'} > /dev/null

salmon_index="$transcripts.salmon"
transcript_name_mapping_file="transcript_name_mapping.txt"

# Ensure that the needed Salmon index exists, exit if not
if [[ ! -d $salmon_index ]]; then
    echo "Expected Salmon index \"$salmon_index\" not found!"
    echo "Please create it first before rerunning this script"
    echo "exiting ..."
    exit 1
fi

# Ensure that the needed transcript name mapping file exists, exit if not
if [[ ! -f $transcript_name_mapping_file ]]; then
    echo "Expected Ensembl to transcript name mapping file \"$transcript_name_mapping_file\" not found!"
    echo "Please create it first before rerunning this script"
    echo "exiting ..."
    exit 1
fi

# Declare array to hold job IDs
declare -a job_ids

function parallel_run_and_wait() {

    ids_file_name="$1"
    command_to_run="$2"

    # Grab the first word of the command to use as the job name
    jobname=$( echo $command_to_run | awk '{print $1}')

    # Run all of the jobs in parallel and capture their job ids in an array
    #   note: the --parsable flag causes sbatch to send bare JOBIDs to standard out (instead of "Submitted batch job JOBID" per JOBID)
    job_ids=($(cat $ids_file_name | parallel "sbatch -o $log_dir/$jobname.%j.out.log --job-name=$jobname --parsable --cpus-per-task=$threads --mem-per-cpu=${mem_per_cpu}G --time=1-00:00:00 --wrap='$command_to_run'" ))

    # Make a comma-delimited list of the previously acquired job ids
    job_id_list=$(IFS=","; echo "${job_ids[*]}")

    # Wait on running jobs before proceeding with the next step
    #   "--wait" flag causes this script to pause until the sbatch job finishes
    sbatch --wait \
           -o $log_dir/wait_job.%j.o.txt  \
           --dependency=$job_id_list \
           --wrap="echo 'Finished waiting for $job_id_list'"
}


# Step 1 - Quality Control

if [ "$run_fastq_qc" == "T" ]; then
    # fastqc runs an analysis for each FASTQ file. Then we run multiqc to aggregate all the FastQC reports
    sbatch --cpus-per-task=$threads \
           -o $log_dir/fastqc.%j.o.txt \
           --mem-per-cpu=${mem_per_cpu}G \
           --wrap="fastqc \
                    --threads=$threads \
                    $fastq_dir/*.f*q.gz \
                    --outdir $qc_dir; \
                    multiqc $qc_dir --outdir $qc_dir"
           #        note: In bash commands, an asterisk matches anything or nothing (but only in the same directory)
           # p.s. No sbatch --wait flag here, since the QC isn't a prerequisite for subsequent steps
fi

# Classify FASTQ reads
#    note: The braces {} will get filled in later with values fed to the "parallel" command

raw_counts_file=$raw_counts_dir/counts.txt
if [ "$reuse_counts" == "F" ]; then
    classification_command=""
    if [ "$is_single" == "T" ]; then

        # Run salmon to classify the reads.
        classification_command="salmon quant \
                                  --libType A `#automatically detect library type` \
                                  --minAssignedFrags $minAssignedFrags \
                                  -i $salmon_index \
                                  -o $salmon_dir/{} \
                                  -r $fastq_dir/{}_R1.f*q.gz"

    else

        # Run salmon to classify the reads.
        classification_command="salmon quant \
                                  --libType A `#automatically detect library type` \
                                  --minAssignedFrags $minAssignedFrags \
                                  -i $salmon_index \
                                  -o $salmon_dir/{} \
                                  -1 $fastq_dir/{}_R1.f*q.gz \
                                  -2 $fastq_dir/{}_R2.f*q.gz"

    fi

    parallel_run_and_wait $ids_file "$classification_command"

    cat $ids_file | python combine.py $salmon_dir | $code_dir/clean_sample_names.sh > $raw_counts_file

    # Step 5 - Check if your samples cluster together with Principal Component Analysis (PCA):
    Rscript $code_dir/pca.r $raw_counts_file > $pca_dir/pca.pdf
fi
