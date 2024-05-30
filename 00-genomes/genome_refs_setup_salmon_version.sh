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

    Usage: $0 --species "<TEXT>" --species "<Species name>" --ref_name "<Reference name>"

    Options:
        --species         In quotes, name of species (defaults to "$species")
        --release         Ensembl release (defaults to $release)
        --ref_name        Reference name (defaults to "$ref_name")

    Example:
        $0 --species "$species" --release $release --ref_name $ref_name

    This will download said Ensembl release, index it with salmon, and create mapping
    of gene IDs to names.

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

# Change species name to all lowercase (following Ensembl conventention for dir names)
species_lc=$( echo $species | tr '[:upper:]' '[:lower:]' )

rm -rf $log_dir
mkdir -p $log_dir

# WARNING: For some less well studied organisms, this is "top_level" instead of "primary_assembly"
genome="$base_ref_name.dna.primary_assembly.fa"

if [[ ! -f $genome ]]; then
    # Download primary assembly (i.e. the one that doesn't include patches and extra haplotypes which would confuse aligners)
    wget http://ftp.ensembl.org/pub/release-$release/fasta/$species_lc/dna/$genome.gz
    pigz -d -c $genome.gz > $genome
fi

if [[ ! -f $annotations ]]; then
    # Get Ensembl annotations
    wget http://ftp.ensembl.org/pub/release-$release/gtf/$species_lc/$annotations.gz
    pigz -d -c $annotations.gz > $annotations
fi

module load salmon     # for indexing transcriptome prior to alignment
module load gffread      # for generating transcriptome from genome and annotation files

transcripts="$base_ref_name.$release.transcripts.fa"
if [[ ! -f $transcripts ]]; then
    # Create transcripts file from genome and annotation files
    # This also indexes the genome
    gffread -w $transcripts -g $genome $annotations
fi

# Get Ensembl transcript ID to gene name mapping file
transcript_name_mapping="$base_ref_name.$release.transcript_name_mapping.txt"
if [[ ! -f $transcript_name_mapping ]]; then
    sbatch --parsable \
        -o $log_dir/gene_names_from_gtf.%j.o.txt \
        --wrap="../00-rcode/transcript_names_from_gtf.sh $annotations > $transcript_name_mapping; ln -s $transcript_name_mapping transcript_name_mapping.txt"
fi

# salmon prep
gentrome="gentrome.$genome.gz"
decoys="decoys.$genome.txt"
if [[ ! -f $gentrome ]]; then
    cat $genome | grep "^>" | cut -d " " -f 1 > $decoys
    sed -i.bak -e 's/>//g' $decoys
    cat $transcripts $genome | gzip -c - > $gentrome
fi

salmon_index="$transcripts.salmon"
mkdir -p $salmon_index

# Index the transcriptome with salmon.
sbatch --parsable \
    -o $log_dir/salmon_index.%j.o.txt \
    --cpus-per-task=12 \
    --mem-per-cpu=6G \
    --wrap="salmon index -t $gentrome -d $decoys --threads 12 -i $salmon_index"

# List loaded modules (and therefore capture their versions)
module list &> module_list.log
