Nucleotide-specific bias correction
-----------------------------------

Reweight high-throughput sequencing reads to account for nucleotide-specific bias
from any source, including assay and sequencing biases.


Overview:

1.  Compute baseline
2.  Compute bias
3.  Correct bias
4.  Plotting and diagnostics


How to use it:

    sh seqbias_pipe.sh [prefix] [bam] [read length] [k] [output]
      prefix        will be appended to all working and result files
      bam           should be aligned and sorted
      read length   reads in BAM file should be the same length
      k             tile size to correct (5 is recommended most of the time)
      output        directory to put all output files
    
    ex. sh seqbias_pipe.sh GM12878Rep1 GM12878Rep1.bam 20 5 results/GM12878Rep1.5mer


Jeremy Wang (Furey Lab)
University of North Carolina at Chapel Hill
November, 2014
