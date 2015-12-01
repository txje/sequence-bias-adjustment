bam=$1
out=$2
chroms=$3

echo "Conversion to wig/bw clips weights to 4x!"
python pipeline/pileup_wig/pileup.py $bam --out $out.wig --clip 4
wigToBigWig $out.wig $chroms $out.bw
