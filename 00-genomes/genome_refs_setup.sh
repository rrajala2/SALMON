#!/usr/bin/bash -l 
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=6G
#SBATCH --time=1-00:00

set -e          # Will EXIT the script if any command returns a non-zero exit status.
set -E          # Make ERR trapping work inside functions too.
set -u          # Variables must be pre-defined before using them.
set -o pipefail # If a pipe fails, returns the error code for the failed pipe even if it isn't the last command in a series of pipes.

# If an error occurs, report the line number to stderr (see https://olivergondza.github.io/2019/10/01/bash-strict-mode.html).
trap 's=$?; echo >&2 "$0: Error on line "$LINENO": $BASH_COMMAND"; exit $s' ERR

# Default reference info
release=107
species="Homo sapiens"
ref_name="GRCh38"

log_dir=logs

function usage() {
    cat <<USAGE

    Usage: $0 --species "<TEXT>" --species "<Specie name>" --ref_name "<Reference name>"

    Options:
        --species         In quotes, name of species (e.g. "Homo sapiens")
        --release         Ensembl release
        --ref_name        Reference name (e.g. "GRCh38" or "GRCm39")

    Example:
        $0 --species "Homo sapiens" --release 106 --ref_name GRCh38 

    This will download said release, index it with hisat2 and kallisto, and also
    create a mapping of gene IDs to names.

USAGE
    exit 1
}

while [ $# -gt 0 ]; do
    case $1 in
        --ref_name)
            ref_name=$2
            ;;
        --release)
            release=$2
            ;;
        --species)
            species=$2
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

species=$( echo "$species" | tr ' ' '_' ) # Replace space with underscore

base_ref_name="$species.$ref_name"
annotations="$base_ref_name.$release.gtf"
species_lc=$( echo $species | tr '[:upper:]' '[:lower:]' )

rm -rf $log_dir
mkdir -p $log_dir

# WARNING: Sometimes this is "top_level" instead of "primary_assembly"
genome="$base_ref_name.dna.primary_assembly.fa"

transcripts="$base_ref_name.$release.transcripts.fa"

# Download primary assembly (i.e. the one that doesn't include patches and extra haplotypes which would confuse aligners)
wget http://ftp.ensembl.org/pub/release-$release/fasta/$species_lc/dna/$genome.gz

# Get Ensembl annotations
wget http://ftp.ensembl.org/pub/release-$release/gtf/$species_lc/$annotations.gz

# Decompress files (for gffread's sake)
pigz -d -c $annotations.gz > $annotations
pigz -d -c $genome.gz > $genome

module load hisat2       # for indexing genome prior to alignment
module load kallisto     # for indexing transcriptome prior to alignment
module load gffread      # for generating transcriptome from genome and annotation files

# Create transcripts file from genome and annotation files
# This also indexes the genome
gffread -w $transcripts -g $genome $annotations

gene_name_mapping="$base_ref_name.$release.gene_name_mapping.txt"
if [[ ! -f $gene_name_mapping ]]; then
    sbatch --parsable \
        -o $log_dir/gene_names_from_gtf.%j.o.txt \
        --wrap="../00-rcode/gene_names_from_gtf.sh $annotations > $gene_name_mapping; ln -s $gene_name_mapping gene_name_mapping.txt"
fi

# Index the transcriptome with kallisto.
sbatch --parsable \
    -o $log_dir/kallisto_index.%j.o.txt \
    --mem=12G \
    --wrap="kallisto index -i $transcripts.kallisto $transcripts"

# Index the reference genome. Needs to be done only once and can be reused.
sbatch --parsable \
    -o $log_dir/hisat2-build.%j.o.txt \
    --mem=12G \
    --cpus-per-task=12 \
    --wrap="hisat2-build --threads=12 $genome $genome"

# Capture variables and create a generic rna_seq script
cat << END > rna_seq_script.sh
#!/bin/bash
genome=$genome
annotations=$annotations
transcripts=$transcripts

level_6_plus_classification.sh \\`echo -e "\n"`
    --fastq-dir fastq \\`echo -e "\n"`
    --ids-file ids.txt \\`echo -e "\n"`
    --control CONTROL \\`echo -e "\n"`
    --condition CONDITION \\`echo -e "\n"`
    --genome 00-genomes/\$genome \\`echo -e "\n"`
    --annotations 00-genomes/\$annotations \\`echo -e "\n"`
    --transcript 00-genomes/\$transcripts \\`echo -e "\n"`

END

# List loaded modules (and therefore capture their versions)
module list &> module_list.log
