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

    Usage: $0 --control <TEXT> --condition <TEXT> --ids-file <FILE> --transcripts <FILE>

    Options:
        --control         Name of control (e.g. "BORED" if control fq files are named BORED_1_R1.fq.gz, BORED_1_R2.fq.gz, BORED_2_R1.fq.gz, BORED_2_R2.fq.gz, BORED_3_R1.fq.gz, BORED_3_R2.fq.gz)
        --condition       Name of condition (e.g. "EXCITED" if control fq files are named EXCITED_1_R1.fq.gz, EXCITED_1_R2.fq.gz, EXCITED_2_R1.fq.gz, EXCITED_2_R2.fq.gz)
        --transcripts     Path to the transcriptome reference FASTA file (e.g. 00-genome/transcripts.fa)
        --ids-file        File containing list of IDs to process

    Example:
        $0 --ids-file ids.txt --control BORED --condition EXCITED --transcripts 00-genome/transcripts.fa

USAGE
    exit 1
}

threads=12    # Number of CPU cores
mem_per_cpu=6 # 6 GBS of RAM per CPU core

# Make some extra directories to help file organization:
code_dir=scripts
raw_counts_dir=04-raw_counts
pca_dir=05-pca
de_dir=06-DE_lists
heatmap_dir=07-heatmaps
log_dir=slurm_logs

mkdir -p $raw_counts_dir
mkdir -p $pca_dir
mkdir -p $de_dir
mkdir -p $heatmap_dir
mkdir -p $log_dir


# Step 0 - getting things set up

while [ $# -gt 0 ]; do
    case $1 in              # check for these flags
        --condition)
            condition=$2    # store value for found flag
            ;;
        --control)
            control=$2
            ;;
        --ids-file)
            ids_file=$2
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
echo ${control?'Please specify the name of the control using --control control_name'} > /dev/null
echo ${condition?'Please specify the name of the condition using --condition condition_name'} > /dev/null
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

# previously created
raw_counts_file=$raw_counts_dir/counts.txt

# Step 6 - Run differential expression analysis with R.
#          Here are the example scripts used by the Biostar handbook for doing DE analysis.
#          Their experimental design is 3 control samples compared to 3 experiment samples. They are
#          testing if there are DE genes between the experiment group and control group.

for de_method in deseq2 edger; do
    base_result_name=${de_method}-${control}_vs_${condition}
    counts_with_ensembl_ids=$de_dir/${base_result_name}.csv
    counts_with_gene_names=$de_dir/${base_result_name}_w_gene_names.tsv

    Rscript $code_dir/$de_method.R $raw_counts_file $control $condition > $counts_with_ensembl_ids

    # We can also generate heatmaps based on the DE gene lists:
    cat $counts_with_ensembl_ids  | Rscript $code_dir/heatmap_from_csv.r > $heatmap_dir/heatmap_${base_result_name}.pdf

    # Get Ensembl ID to gene name mapping file
    $code_dir/add_gene_names.sh $counts_with_ensembl_ids $transcript_name_mapping_file >  $counts_with_gene_names

    cat $counts_with_gene_names | Rscript $code_dir/HEATMAP_new.R > $heatmap_dir/heatmap_${base_result_name}_w_gene_names.pdf
done
