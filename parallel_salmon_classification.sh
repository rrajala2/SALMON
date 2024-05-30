#!/usr/bin/bash -l 
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=6G
#SBATCH --time=1-00:00
#SBATCH --mail-user=$USER@omrf.org
#SBATCH --mail-type=END,FAIL

config_file=config.csv

transcripts_file=$1

./salmon_align.sh --transcripts $transcripts_file --fastq-dir 00-reads --ids-file ids.txt

job_id=$( sbatch --parsable --wrap="cat ids.txt | scripts/combine_tpm.py 02-salmon > 02-salmon/tpm_counts.tsv")

sbatch --dependency=$job_id --wrap="scripts/add_gene_names.sh 02-salmon/tpm_counts.tsv transcript_name_mapping.txt > 02-salmon/tpm_counts_w_gene_names.tsv"

cat $config_file | parallel --colsep=, "./salmon_de.sh --transcripts $transcripts_file --ids-file ids.txt --control {1} --condition {2}"
