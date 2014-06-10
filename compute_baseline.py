# ignore alignability for the moment, we'll assume everything within 100bp of a read is alignable

import numpy
import argparse
import pysam

SAMPLE_MARGIN = 25 # on either side
INVERSE = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
CHROMS = ['chr%i' % i for i in xrange(1,23)] + ['chrX', 'chrY']

def readWig(wig_file):
  fin = open(wig_file, 'r')
  masks = {}
  mask = []
  chrom = None
  start = 1
  count = 0
  for line in fin:
    if line[0] == 'f': # header line ex. "fixedStep chrom=chr1 step=1 ..."
      params = {param.split('=')[0]:param.split('=')[1] for param in line.strip().split(' ') if '=' in param}
      if params["chrom"] != chrom:
        if len(mask) > 0:
          masks[chrom] = numpy.array(mask, dtype='b') # boolean?
          mask = []
        chrom = params["chrom"]
        print chrom
        '''
        if chrom == "chr10":
          break
        '''
      elif start + count != int(params["start"]):
        print "GAP!"
      start = int(params["start"])
      count = 0
    else:
      mask.append(int(float(line.strip())))
      count += 1
  if len(mask) > 0:
    masks[chrom] = numpy.array(mask, dtype='b') # boolean?
    mask = []
  fin.close()

  # fill out the empty chroms for testing...
  for c in xrange(len(CHROMS)):
    if not masks.has_key(CHROMS[c]):
      print "missing", CHROMS[c]
      masks[c] = numpy.zeros(0, dtype='b')
    else:
      masks[c] = masks[CHROMS[c]]
      del masks[CHROMS[c]]

  return masks

def getRandom(chrom, first, last, cnt, alignability_mask):
  if numpy.sum(alignability_mask[chrom][first:last]) > 0:
    nums = numpy.random.choice([i for i in xrange(first, last) if alignability_mask[chrom][i] == 1], cnt, replace=True)
  else:
    nums = numpy.random.randint(first, last, cnt)
  if not isinstance(nums, list):
    nums = [nums]
  return nums

def main(bam_npy_file, fasta_file, k, output_file, exp=1, limit=None, window_max=SAMPLE_MARGIN*2):

  print "Reading alignability mask..."
  alignability_mask = readWig("filter/wgEncodeDukeMapabilityUniqueness20bp.wig")

  baseline_kmer_counts = {}

  trailing = 0
  leading = 1
  density_count = 1

  print "Reading FASTA seqs..."
  fa = pysam.Fastafile(fasta_file)
  refs = [fa.fetch(c).upper() for c in CHROMS]
  ref_lens = [len(r) for r in refs]

  print "Reading BAM..."
  bam = numpy.load(bam_npy_file)
  tot_reads = len(bam)
  print "%i reads" % tot_reads
  num_reads = 0
  b = 0
  while b < tot_reads:
    read = bam[b]

    if b % 10**6 == 0:
      print "%i (%.2f%%)" % (b, float(b)/tot_reads*100)

    while leading < tot_reads and (bam[leading][0] < read[0] or (bam[leading][0] == read[0] and bam[leading][1] <= read[1] + SAMPLE_MARGIN)):
      if bam[leading][1] != bam[leading-1][1]: # count only 1 read at each unique position
        density_count += 1
      leading += 1
    while trailing < b and (bam[trailing][0] < read[0] or bam[trailing][1] < read[1] - SAMPLE_MARGIN):
      if bam[trailing][1] != bam[trailing-1][1]:
        density_count -= 1
      trailing += 1

    if density_count < 0:
      print density_count, "too low"

    first = max(0, read[1] - SAMPLE_MARGIN) # inclusive
    last = min(ref_lens[read[0]] - k, read[1] + SAMPLE_MARGIN + 1 - k) # exclusive
    for pos in getRandom(read[0], first, last, min(int(density_count ** exp), window_max), alignability_mask): # get (density ** exp) random positions

      # get sequence
      if read[2]:
        seq = ''.join([INVERSE[a] for a in refs[read[0]][pos : pos + k][::-1]])
      else:
        seq = refs[read[0]][pos : pos + k]

      # count k-mers
      #for i in xrange(len(seq) - k + 1):
      baseline_kmer_counts[seq] = baseline_kmer_counts.get(seq, 0) + 1

    b = leading + 1

  fout = open(output_file, 'w')
  kmers = baseline_kmer_counts.keys()
  total_kmers = sum(v for k,v in baseline_kmer_counts.iteritems())
  fout.write(','.join(kmers) + "\n")
  fout.write(','.join(['%.4f' % (float(baseline_kmer_counts[k])/total_kmers) for k in kmers]))
  fout.close()

if __name__ == "__main__":
  parser = argparse.ArgumentParser(description = "Compute baseline k-mer frequencies")
  parser.add_argument("bam", help="BAM.npy file")
  parser.add_argument("ref", help="Fasta file")
  parser.add_argument("k", help="k-mer size", type=int)
  parser.add_argument("out", help="Output (csv) file")
  parser.add_argument("--exp", help="Density exponent", type=float, default=1)
  parser.add_argument("--limit", help="Reads to look at, tiled across entire BAM", type=int)
  parser.add_argument("--window_max", help="Maximum # reads to count in a single window", type=int)
  args = parser.parse_args()

  main(args.bam, args.ref, args.k, args.out, args.exp, args.limit, args.window_max)
