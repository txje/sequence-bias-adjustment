Correcting nucleotide-specific bias in high-throughput sequencing data
----------------------------------------------------------------------

Reweight high-throughput sequencing reads to account for nucleotide-specific bias
from any source, including assay and sequencing biases.


Overview:

1.  Compute baseline
2.  Compute bias
3.  Correct bias
4.  Plotting and diagnostics


How to use it:

    sh seqbias_pipe.sh <ref> <chroms> <prefix> <bam> <k> <outdir> [--resume <step>]
      ref           reference genome FASTA file
      path          prefix of data files
      chroms        list of chromosomes (one per line) to correct
      prefix        will be appended to all working and result files
      bam           should be aligned and sorted
      k             tile size to correct (5 is recommended most of the time)
      outdir        directory to put all output files


Jeremy Wang, Ph.D.
Department of Genetics
University of North Carolina at Chapel Hill

MIT License (see LICENSE.txt)
