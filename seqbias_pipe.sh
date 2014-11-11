timestamp() {
  date +"%Y%m%d%H%M%S"
}

prefix=$1
bam=$2
read_len=$3
k=$4
outdir=$5
bamfiltered=data/$prefix.filtered.bam
bamsorted=data/$prefix.filtered.sorted.bam
bamnpy=$bamfiltered.npy
ref=/proj/fureylab/genomes/human/hg19_reference/hg19/hg19.fa
baseline=$outdir/$prefix.read_50_baseline.csv
kmer_bias=$outdir/$prefix.${k}mer_frequencies.csv
tile_cov=$outdir/$prefix.tile_covariance.npy
corrected_weights=$outdir/$prefix.${k}mer_adjusted.read_weights.npy

set -e

mkdir -p data
mkdir --parents $outdir

# ---------------------------------------------------------
# filter, index, and convert BAM to more usable format
#
if [ ! -e $bamnpy ]
then
  echo [$(timestamp)] start: BAM filter, index, convert to npy
  sh filter/filterBlacklistBam.sh $bam $bamfiltered
  samtools sort $bamfiltered data/$prefix.filtered.sorted # will append the .bam
  mv $bamsorted $bamfiltered
  samtools index $bamfiltered
  python bam2npy.py $bamfiltered
  echo [$(timestamp)] end: BAM filter, index, convert to npy
fi
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
# compute k-mer bias
#
echo [$(timestamp)] start: Compute $k-mer bias
python compute_bias.py $bamnpy $ref $k $kmer_bias --read_len $read_len
echo [$(timestamp)] end: Compute $k-mer bias
# ---------------------------------------------------------

# ---------------------------------------------------------
# compute tile covariance matrix
#
echo [$(timestamp)] start: Compute tile covariance matrix
python correlate_bias.py $bamnpy $ref $kmer_bias $tile_cov
echo [$(timestamp)] end: Compute tile covariance matrix
# ---------------------------------------------------------

# ---------------------------------------------------------
# correct using k-mer bias and tile covariance matrix
#
echo [$(timestamp)] start: Correct bias
python correct_bias.py $bamnpy $ref $baseline $kmer_bias $outdir/$prefix.${k}mer_adjusted.allele_frequencies.csv $corrected_weights $tile_cov --read_len $read_len
echo [$(timestamp)] end: Correct bias
# ---------------------------------------------------------

# ---------------------------------------------------------
# compute corrected k-mer bias
#
echo [$(timestamp)] start: Compute corrected $k-mer bias
python compute_bias.py $corrected_weights $ref $k $outdir/$prefix.${k}mer_adjusted.${k}mer_frequencies.csv --read_len $read_len
echo [$(timestamp)] end: Compute corrected $k-mer bias
# ---------------------------------------------------------

# ---------------------------------------------------------
# various plots, results, diagnostics
#
echo [$(timestamp)] start: Plot adjusted nucleotide frequencies
python plot_csv.py $outdir/$prefix.allele_frequencies.csv --out $outdir/$prefix.allele_frequencies.png
python plot_csv.py $outdir/$prefix.${k}mer_adjusted.allele_frequencies.csv --out $outdir/$prefix.${k}mer_adjusted.allele_frequencies.png
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
rm -f last
dir=$(ls -t | head -n 1)
ln -s $dir last
