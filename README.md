Correcting nucleotide-specific bias in high-throughput sequencing data
----------------------------------------------------------------------

Reweight high-throughput sequencing reads to account for nucleotide-specific bias
from any source, including assay and sequencing biases.

Requirements:
* samtools (in PATH)
* Python (2.7)
  * numpy
  * pysam
  * matplotlib (pyplot)

Installation:

    git clone http://github.com/txje/sequence-bias-adjustment
    cd sequence-bias-adjustment

Usage:

    sh seqbias_pipe.sh <ref> <chroms> <bam> <k> <outdir> <prefix> [--resume <step>]
      ref           reference genome FASTA file
      chroms        list of chromosomes (one per line) to correct
      bam           should be aligned and blacklist filtered (if needed)
      k             tile size to correct (5 is recommended most of the time)
      outdir        directory to put all output files
      prefix        prepended to all intermediate and output files

Example -- to run bias correction on an example ENCODE DNase-seq data set:

    git clone http://github.com/txje/sequence-bias-adjustment
    cd sequence-bias-adjustment
    mkdir example_data
    cd example_data
    wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeOpenChromDnase/wgEncodeOpenChromDnaseGm12878AlnRep1.bam
    wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeOpenChromDnase/wgEncodeOpenChromDnaseGm12878AlnRep1.bam.bai
    wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes
    wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit
    wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa
    chmod 700 twoBitToFa
    ./twoBitToFa hg19.2bit hg19.fa
    cd ..
    sh seqbias_pipe.sh example_data/hg19.fa example_data/hg19.chrom.sizes example_data/wgEncodeOpenChromDnaseGm12878AlnRep1.bam 5 example_results gm12878.dnase


Jeremy Wang, Ph.D.  
Department of Genetics  
University of North Carolina at Chapel Hill

MIT License (see LICENSE.txt)
