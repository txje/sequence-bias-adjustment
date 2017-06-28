timestamp() {
  date +"%Y%m%d%H%M%S"
}

ref=$1
chroms=$2
bam=$3
k=$4
outdir=$5
prefix=$6
opt=$7
resume=0
if [[ $opt == "--resume" ]]
then
  step=$8
  case $step in
    baseline)
      resume=1;;
    bias)
      resume=2;;
    covariance)
      resume=3;;
    correct)
      resume=4;;
    rebias)
      resume=5;;
    bam)
      resume=6;;
    bigwig)
      resume=7;;
  esac
  echo "Resuming at step $resume ($step)"
fi
bamnpy=$outdir/$prefix.bam.npy

baseline=$outdir/$prefix.read_50_baseline.csv
kmer_bias=$outdir/$prefix.${k}mer_frequencies.csv
tile_cov=$outdir/$prefix.tile_covariance.npy
corrected_weights=$outdir/$prefix.${k}mer_adjusted.read_weights.npy

set -e

mkdir -p $outdir

# ---------------------------------------------------------
# filter, index, and convert BAM to more usable format
#
if [ ! -e $bamnpy ]
then
  echo [$(timestamp)] start: BAM filter, index, convert to npy
  echo step begin

  read_len=$(samtools view $bam | head -n 1 | awk '{print $10}' - | wc | awk '{print $3}' -)

  python bam2npy.py $bam $chroms $bamnpy
  echo [$(timestamp)] end: BAM filter, index, convert to npy

else
  read_len=$(samtools view $bam | head -n 1 | awk '{print $10}' - | wc | awk '{print $3}' -)
fi

# get read length
echo "$prefix read length: $read_len"
if [ $read_len == 2 ]
then
  read_len=20
  echo "$prefix read length adjusted to: $read_len"
fi

echo "$prefix read length: $read_len" >> debug.out
# ---------------------------------------------------------

# ---------------------------------------------------------
# compute nucleotide bias
#
if [ $resume -lt 1 ]
then
  echo [$(timestamp)] start: Compute nucleotide bias
  python compute_bias.py $bamnpy $ref $chroms 1 $outdir/$prefix.allele_frequencies.csv --read_len $read_len
  echo [$(timestamp)] end: Compute nucleotide bias
fi
# ---------------------------------------------------------

# ---------------------------------------------------------
# compute baseline (read +- 50)
#
if [ $resume -lt 2 ]
then
  echo [$(timestamp)] start: Compute $k-mer baseline
  echo step baseline
  python compute_baseline.py $bamnpy $ref $chroms $k $baseline
  # --mask $alignability_mask
  echo [$(timestamp)] end: Compute $k-mer baseline
fi
# ---------------------------------------------------------

# ---------------------------------------------------------
# compute k-mer bias
#
if [ $resume -lt 3 ]
then
  echo [$(timestamp)] start: Compute $k-mer bias
  echo step bias
  python compute_bias.py $bamnpy $ref $chroms $k $kmer_bias --read_len $read_len
  echo [$(timestamp)] end: Compute $k-mer bias
fi
# ---------------------------------------------------------

# ---------------------------------------------------------
# compute tile covariance matrix
#
if [ $resume -lt 4 ]
then
  echo [$(timestamp)] start: Compute tile covariance matrix
  echo step covariance
  python correlate_bias.py $bamnpy $ref $chroms $kmer_bias $tile_cov --plot
  echo [$(timestamp)] end: Compute tile covariance matrix
fi
# ---------------------------------------------------------

# ---------------------------------------------------------
# correct using k-mer bias and tile covariance matrix
#
if [ $resume -lt 5 ]
then
  echo [$(timestamp)] start: Correct bias
  echo step correct
  python correct_bias.py $bamnpy $ref $chroms $baseline $kmer_bias $outdir/$prefix.${k}mer_adjusted.allele_frequencies.csv $corrected_weights $tile_cov --read_len $read_len
  echo [$(timestamp)] end: Correct bias
fi
# ---------------------------------------------------------

if [ $resume -lt 6 ]
then
# ---------------------------------------------------------
# compute corrected nucleotide bias
#
  echo [$(timestamp)] start: Compute corrected nucleotide bias
  echo step rebias
  python compute_bias.py $corrected_weights $ref $chroms 1 $outdir/$prefix.${k}mer_adjusted.allele_frequencies.csv --read_len $read_len
  echo [$(timestamp)] end: Compute corrected nucleotide bias
# ---------------------------------------------------------

# ---------------------------------------------------------
# compute corrected k-mer bias
#
  echo [$(timestamp)] start: Compute corrected $k-mer bias
  python compute_bias.py $corrected_weights $ref $chroms $k $outdir/$prefix.${k}mer_adjusted.${k}mer_frequencies.csv --read_len $read_len
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
fi

# ---------------------------------------------------------
# Output BAM format with XW weight flag
#
if [ $resume -lt 7 ]
then
  echo [$(timestamp)] start: Writing BAM output
  echo step bam
  python npy2bam.py $chroms $outdir/$prefix.${k}mer_adjusted.read_weights.npy $bam $outdir/$prefix.adjusted.bam --noy --tag
  echo [$(timestamp)] end: Writing BAM output

  echo [$(timestamp)] start: Indexing BAM
  samtools index $outdir/$prefix.adjusted.bam
  echo [$(timestamp)] end: Indexing BAM
fi
# ---------------------------------------------------------

# ---------------------------------------------------------
# Convert weighted bam to wig/bw
#
#if [ $resume -lt 8 ]
#then
#  echo [$(timestamp)] start: Writing wig/bw output
#  echo step bigwig
#  # right now, this is fixed to hg19
#  sh pileup_wig/wigify.sh $outdir/$prefix.adjusted.bam $outdir/$prefix.adjusted $chroms
#  echo [$(timestamp)] end: Writing wig/bw output
#fi
# ---------------------------------------------------------

cd $outdir
cd ..
rm -f last
dir=$(ls -t | head -n 1)
ln -s $dir last
