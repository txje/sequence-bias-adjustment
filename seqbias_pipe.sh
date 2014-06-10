timestamp() {
  date +"%Y%m%d%H%M%S"
}

prefix=$1
bam=$2
read_len=$3
k=$4
outdir=$5
bamnpy=$prefix.filtered.bam.npy
ref=/proj/fureylab/genomes/human/hg19_reference/hg19/hg19.fa
baseline=$outdir/$prefix.read_50_baseline.csv
# random baseline (suboptimal)
#baseline=hg19_baseline_alignable_5mer_frequency.csv
kmer_bias=$outdir/$prefix.${k}mer_frequencies.csv

set -e

mkdir --parents $outdir

# ---------------------------------------------------------
# filter, index, and convert BAM to more usable format
#
#echo [$(timestamp)] start: BAM filter, index, convert to npy
#sh filter/filterBlacklistBam.sh $bam $prefix.filtered.bam
#samtools index $prefix.filtered.bam
#python filter/bam2npy.py $prefix.filtered.bam
#echo [$(timestamp)] end: BAM filter, index, convert to npy
# ---------------------------------------------------------

# ---------------------------------------------------------
# compute nucleotide bias
#
echo [$(timestamp)] start: Compute nucleotide bias
python compute_bias.py $bamnpy $ref 1 $outdir/$prefix.allele_frequencies.csv --read_len $read_len
echo [$(timestamp)] end: Compute nucleotide bias
# ---------------------------------------------------------

# ---------------------------------------------------------
# compute baseline (read +- 50)
#
echo [$(timestamp)] start: Compute $k-mer baseline
python compute_baseline.py $bamnpy $ref $k $baseline
echo [$(timestamp)] end: Compute $k-mer baseline
# ---------------------------------------------------------

# ---------------------------------------------------------
# compute and correct for bias using k-mer method
#
echo [$(timestamp)] start: Compute $k-mer bias
python compute_bias.py $bamnpy $ref $k $kmer_bias --read_len $read_len
echo [$(timestamp)] end: Compute $k-mer bias

echo [$(timestamp)] start: Correct bias
python correct_bias.py $bamnpy $ref $baseline $kmer_bias $outdir/$prefix.${k}mer_adjusted.allele_frequencies.csv $outdir/$prefix.${k}mer_adjusted.read_weights.csv.npy --read_len $read_len
echo [$(timestamp)] end: Correct bias
# ---------------------------------------------------------

# ---------------------------------------------------------
# various plots, results, diagnostics
#
echo [$(timestamp)] start: Plot adjusted nucleotide frequencies
python plot_csv.py $outdir/$prefix.allele_frequencies.csv $outdir/$prefix.allele_frequencies.png
python plot_csv.py $outdir/$prefix.${k}mer_adjusted.allele_frequencies.csv $outdir/$prefix.${k}mer_adjusted.allele_frequencies.png
echo [$(timestamp)] end: Plot adjusted nucleotide frequencies
#
# elsewhere:
# - Compute and plot TFBS bias (before and after correction)
#    see bsub_tfbs.sh
# - Correlation with ChIP-seq peaks
#    see bsub_compare.sh
# ---------------------------------------------------------

cd $outdir
cd ..
dir=$(ls -t | head -n 1)
ln -s $dir last
